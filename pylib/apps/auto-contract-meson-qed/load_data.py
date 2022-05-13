import qlat as q
import rbc_ukqcd_params as rup
import rbc_ukqcd as ru

from auto_contractor.eval import *

import os

from jobs import *

def get_prop_wsrc(prop_cache, inv_type, t_src, tag_snk_type):
    cache_type_dict = {
            "wsrc_wsnk ; psel_ts": "psel_ts",
            "wsrc ; fsel": "fsel",
            "wsrc ; psel": "psel",
            }
    cache_type = cache_type_dict[tag_snk_type]
    prop_cache_prob = prop_cache["prob"]
    prop_cache_type = prop_cache[cache_type]
    tslice = t_src
    def mk_tag(inv_acc):
        return f"tslice={tslice} ; type={inv_type} ; accuracy={inv_acc} ; {tag_snk_type}"
    tag = mk_tag(inv_acc = 1)
    tag1 = mk_tag(inv_acc = 2)
    prob = prop_cache_prob[f"type={inv_type} ; accuracy=1 ; wsrc ; prob"]
    # level light_accuracy strange_accuracy
    # 0     inv_acc=1      inv_acc=1
    # 3     inv_acc=2      inv_acc=2
    assert prob == 1
    val = prop_cache_type.get(tag)
    if tag1 not in prop_cache_type:
        return val
    source_specification = repr(("wall", t_src,))
    val_list = [ val, prop_cache_type.get(tag1), ]
    rel_acc_list = [ 0, 3, ]
    prob_list = [ 1, prop_cache_prob[f"type={inv_type} ; accuracy=2 ; wsrc ; prob"], ]
    return mk_ama_val(val, source_specification, val_list, rel_acc_list, prob_list)

def get_prop_wsnk_wsrc(prop_cache, inv_type, t_snk, t_src):
    sp_prop = get_prop_wsrc(prop_cache, inv_type, t_src, "wsrc_wsnk ; psel_ts")
    def f(x):
        return x.get_elem(t_snk)
    return ama_apply1(f, sp_prop)

def get_prop_psnk_wsrc_fsel(prop_cache, inv_type, xg_snk, t_src, fselc_pos_dict):
    assert isinstance(xg_snk, tuple) and len(xg_snk) == 4
    idx_snk = fselc_pos_dict[xg_snk]
    def f(x):
        return x.get_elem(idx_snk)
    return ama_apply1(f, get_prop_wsrc(prop_cache, inv_type, t_src, "wsrc ; fsel"))

def get_prop_psnk_wsrc_psel(prop_cache, inv_type, xg_snk, t_src, psel_pos_dict):
    assert isinstance(xg_snk, tuple) and len(xg_snk) == 4
    idx_snk = psel_pos_dict[xg_snk]
    def f(x):
        return x.get_elem(idx_snk)
    return ama_apply1(f, get_prop_wsrc(prop_cache, inv_type, t_src, "wsrc ; psel"))

### -------

def get_prop_psrc(prop_cache, inv_type, xg_src, tag_snk_type):
    cache_type_dict = {
            "psrc_wsnk ; psel_ts": "psel_ts",
            "psrc ; fsel": "fsel",
            "psrc ; psel": "psel",
            }
    cache_type = cache_type_dict[tag_snk_type]
    prop_cache_prob = prop_cache["prob"]
    prop_cache_type = prop_cache[cache_type]
    xg = xg_src
    xg_str = f"({xg[0]},{xg[1]},{xg[2]},{xg[3]})"
    def mk_tag(inv_acc):
        return f"xg={xg_str} ; type={inv_type} ; accuracy={inv_acc} ; {tag_snk_type}"
    tag = mk_tag(inv_acc = 0)
    tag1 = mk_tag(inv_acc = 1)
    tag2 = mk_tag(inv_acc = 2)
    prob = prop_cache_prob[f"type={inv_type} ; accuracy=0 ; psrc ; prob"]
    # level light_accuracy strange_accuracy
    # 0     inv_acc=0      inv_acc=0/zero_prop
    # 1     inv_acc=0      inv_acc=0
    # 2     inv_acc=1      inv_acc=1
    # 3     inv_acc=2      inv_acc=2
    source_specification = repr(("point", tuple(xg_src),))
    if prob == 1.0:
        val = prop_cache_type.get(tag)
        assert val is not None
        if tag1 not in prop_cache_type:
            return val
    else:
        assert prob < 1
        assert inv_type == 1
        val = 0
        if tag not in prop_cache_type:
            return val
    val_list = [
            val,
            prop_cache_type.get(tag),
            prop_cache_type.get(tag1),
            prop_cache_type.get(tag2),
            ]
    rel_acc_list = [ 0, 1, 2, 3, ]
    # Should use the same prob list for strange and light for correct AMA merge!!!
    # It should be the same for accuracy=1,2
    # At present, it is a hack for the accuracy=0 case, where number strange quark props may be less than the number of light quark props.
    prob_list = [
            1.0,
            prop_cache_prob[f"type=1 ; accuracy=0 ; psrc ; prob"],
            prop_cache_prob[f"type=1 ; accuracy=1 ; psrc ; prob"],
            prop_cache_prob[f"type=1 ; accuracy=2 ; psrc ; prob"],
            ]
    return mk_ama_val(val, source_specification, val_list, rel_acc_list, prob_list)

def get_prop_wsnk_psrc(prop_cache, inv_type, t_snk, xg_src):
    assert isinstance(xg_src, tuple) and len(xg_src) == 4
    sp_prop = get_prop_psrc(prop_cache, inv_type, xg_src, "psrc_wsnk ; psel_ts")
    def f(x):
        if isinstance(x, int) and x == 0:
            return 0
        return x.get_elem(t_snk)
    return ama_apply1(f, sp_prop)

def get_prop_psnk_psrc_fsel(prop_cache, inv_type, xg_snk, xg_src, fselc_pos_dict):
    assert isinstance(xg_src, tuple) and len(xg_src) == 4
    assert isinstance(xg_snk, tuple) and len(xg_snk) == 4
    idx_snk = fselc_pos_dict[xg_snk]
    def f(x):
        if isinstance(x, int) and x == 0:
            return 0
        return x.get_elem(idx_snk)
    return ama_apply1(f, get_prop_psrc(prop_cache, inv_type, xg_src, "psrc ; fsel"))

def get_prop_psnk_psrc_psel(prop_cache, inv_type, xg_snk, xg_src, psel_pos_dict):
    assert isinstance(xg_src, tuple) and len(xg_src) == 4
    assert isinstance(xg_snk, tuple) and len(xg_snk) == 4
    idx_snk = psel_pos_dict[xg_snk]
    def f(x):
        if isinstance(x, int) and x == 0:
            return 0
        return x.get_elem(idx_snk)
    return ama_apply1(f, get_prop_psrc(prop_cache, inv_type, xg_src, "psrc ; psel"))

### -------

def get_prop_rand_u1_fsel(prop_cache, inv_type):
    inv_acc = 2
    tag = f"type={inv_type} ; accuracy={inv_acc} ; rand_u1 ; fsel"
    return prop_cache["fsel"].get(tag)

def get_prop_psnk_rand_u1_fsel(prop_cache, inv_type, xg_snk, fsel_pos_dict):
    assert isinstance(xg_snk, tuple) and len(xg_snk) == 4
    idx_snk = fsel_pos_dict[xg_snk]
    def f(x):
        return x.get_elem(idx_snk)
    return ama_apply1(f, get_prop_rand_u1_fsel(prop_cache, inv_type))

### -------

@q.timer
def get_prop_snk_src(prop_cache, flavor, p_snk, p_src, *, psel_pos_dict, fsel_pos_dict, fselc_pos_dict):
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
        msc = ama_apply1(as_mspincolor, msc)
    elif type_snk[:5] == "point" and type_src == "wall":
        assert isinstance(pos_src, int)
        if type_snk == "point":
            assert pos_snk in psel_pos_dict
            msc = get_prop_psnk_wsrc_psel(prop_cache, inv_type, pos_snk, pos_src, psel_pos_dict)
        else:
            assert pos_snk in fsel_pos_dict
            msc = get_prop_psnk_wsrc_fsel(prop_cache, inv_type, pos_snk, pos_src, fselc_pos_dict)
        msc = ama_apply1(as_mspincolor, msc)
    elif type_snk == "wall" and type_src[:5] == "point":
        assert isinstance(pos_snk, int)
        if type_src == "point":
            assert pos_src in psel_pos_dict
            msc = get_prop_psnk_wsrc_psel(prop_cache, inv_type, pos_src, pos_snk, psel_pos_dict)
        else:
            assert pos_src in fsel_pos_dict
            msc = get_prop_psnk_wsrc_fsel(prop_cache, inv_type, pos_src, pos_snk, fselc_pos_dict)
        msc = ama_apply1(as_mspincolor, msc)
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
            msc = get_prop_psnk_psrc_psel(prop_cache, inv_type, pos_snk, pos_src, psel_pos_dict)
            msc = ama_apply1(as_mspincolor, msc)
        elif type_src == "point":
            # means we use point source at the source location
            assert pos_src in psel_pos_dict
            assert pos_snk in fsel_pos_dict
            msc = get_prop_psnk_psrc_fsel(prop_cache, inv_type, pos_snk, pos_src, fselc_pos_dict)
            msc = ama_apply1(as_mspincolor, msc)
        elif type_snk == "point":
            # means we use point source at the sink location
            assert pos_snk in psel_pos_dict
            assert pos_src in fsel_pos_dict
            msc = get_prop_psnk_psrc_fsel(prop_cache, inv_type, pos_src, pos_snk, fselc_pos_dict)
            msc = ama_apply1(as_mspincolor, msc)
            msc = ama_apply1(g5_herm, msc)
        elif pos_snk == pos_src and flavor in rand_u1_flavors:
            # use the rand_u1 source
            assert pos_snk in fsel_pos_dict
            msc = get_prop_psnk_rand_u1_fsel(prop_cache, inv_type, pos_snk, fsel_pos_dict)
            msc = ama_apply1(as_mspincolor, msc)
        else:
            # if nothing else work, try use point src propagator
            assert False
    else:
        raise Exception("get_prop_snk_src unknown p_snk={p_snk} p_src={p_src}")
    return msc

### -------

@q.timer
def load_prop_wsrc_psel(job_tag, traj, flavor, *, wi, psel, fsel, fselc, gt):
    # cache_psel[f"tslice={tslice} ; type={inv_type} ; accuracy={inv_acc} ; wsrc ; psel"]
    # cache_psel_ts[f"tslice={tslice} ; type={inv_type} ; accuracy={inv_acc} ; wsrc_wsnk ; psel_ts"]
    # cache_prob[f"type={inv_type} ; accuracy={inv_acc} ; wsrc ; prob"]
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
    path_sp = f"psel-prop-wsrc-{flavor_tag}/{job_tag}/traj={traj}"
    gt_inv = gt.inv()
    count = { 1: 0, 2: 0, }
    for idx, tslice, inv_type, inv_acc in wi:
        if inv_type != flavor_inv_type:
            continue
        q.displayln_info(0, f"load_prop_wsrc_psel: idx={idx} tslice={tslice} inv_type={inv_type} path_sp={path_sp}")
        tag = f"tslice={tslice} ; type={inv_type} ; accuracy={inv_acc}"
        # load psel psnk prop
        fn_sp = os.path.join(path_sp, f"{tag}.lat")
        sp_prop = q.PselProp(psel)
        sp_prop.load(get_load_path(fn_sp))
        sp_prop = gt_inv * sp_prop
        cache_psel[f"{tag} ; wsrc ; psel"] = sp_prop
        # load wsnk prop
        fn_spw = os.path.join(path_sp, f"{tag} ; wsnk.lat")
        spw_prop = q.PselProp(psel_ts)
        spw_prop.load(get_load_path(fn_spw))
        cache_psel_ts[f"{tag} ; wsrc_wsnk ; psel_ts"] = spw_prop
        count[inv_acc] += 1
    assert count[1] == total_site[3]
    cache_prob[f"type={flavor_inv_type} ; accuracy=1 ; wsrc ; prob"] = 1
    cache_prob[f"type={flavor_inv_type} ; accuracy=2 ; wsrc ; prob"] = get_prob_exact_wsrc(job_tag)

@q.timer
def load_prop_wsrc_fsel(job_tag, traj, flavor, *, wi, psel, fsel, fselc, gt):
    # need to load psel first
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
    gt_inv = gt.inv()
    count = { 1: 0, 2: 0, }
    for idx, tslice, inv_type, inv_acc in wi:
        if inv_type != flavor_inv_type:
            continue
        q.displayln_info(0, f"load_prop_wsrc_fsel: idx={idx} tslice={tslice} inv_type={inv_type} path_s={path_s}")
        tag = f"tslice={tslice} ; type={inv_type} ; accuracy={inv_acc}"
        # load fsel psnk prop
        sc_prop = q.SelProp(fselc)
        sc_prop.load_double_from_float(sfr, tag)
        sc_prop = gt_inv * sc_prop
        cache_fsel[f"{tag} ; wsrc ; fsel"] = sc_prop
        # check psel psnk prop
        sp_prop_diff = q.PselProp(psel)
        sp_prop_diff @= cache_fsel[f"{tag} ; wsrc ; fsel"]
        sp_prop_diff -= cache_psel[f"{tag} ; wsrc ; psel"]
        assert sp_prop_diff.qnorm() <= 1e-14 * cache_psel[f"{tag} ; wsrc ; psel"].qnorm()
        # increase count
        count[inv_acc] += 1
    sfr.close()
    assert count[1] == total_site[3]
    assert cache_prob[f"type={flavor_inv_type} ; accuracy=1 ; wsrc ; prob"] == 1
    assert cache_prob[f"type={flavor_inv_type} ; accuracy=2 ; wsrc ; prob"] == get_prob_exact_wsrc(job_tag)

@q.timer
def load_prop_psrc_psel(job_tag, traj, flavor, *, psel, fsel, fselc):
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
    path_sp = f"psel-prop-psrc-{flavor_tag}/{job_tag}/traj={traj}"
    count = { 0: 0, 1: 0, 2: 0, }
    inv_type = flavor_inv_type
    idx = 0
    xg_list = psel.to_list()
    for xg, inv_acc in [ (xg, inv_acc) for xg in xg_list for inv_acc in (0, 1, 2,) ]:
        xg_str = f"({xg[0]},{xg[1]},{xg[2]},{xg[3]})"
        tag = f"xg={xg_str} ; type={inv_type} ; accuracy={inv_acc}"
        fn_sp = os.path.join(path_sp, f"{tag}.lat")
        fn_sp_load = get_load_path(fn_sp)
        if fn_sp_load is None:
            continue
        q.displayln_info(0, f"load_prop_psrc_psel: idx={idx} ; {tag} ; path_sp={path_sp}")
        idx += 1
        # load psel psnk prop
        sp_prop = q.PselProp(psel)
        sp_prop.load(fn_sp_load)
        cache_psel[f"{tag} ; psrc ; psel"] = sp_prop
        # load wsnk prop
        fn_spw = os.path.join(path_sp, f"{tag} ; wsnk.lat")
        fn_spw_load = get_load_path(fn_spw)
        if fn_spw_load is not None:
            spw_prop = q.PselProp(psel_ts)
            spw_prop.load(fn_spw_load)
            cache_psel_ts[f"{tag} ; psrc_wsnk ; psel_ts"] = spw_prop
        count[inv_acc] += 1
    cache_prob[f"type={flavor_inv_type} ; accuracy=0 ; psrc ; prob"] = count[0] / len(xg_list)
    cache_prob[f"type={flavor_inv_type} ; accuracy=1 ; psrc ; prob"] = rup.dict_params[job_tag]["prob_acc_1_psrc"]
    cache_prob[f"type={flavor_inv_type} ; accuracy=2 ; psrc ; prob"] = rup.dict_params[job_tag]["prob_acc_2_psrc"]

@q.timer
def load_prop_psrc_fsel(job_tag, traj, flavor, *, psel, fsel, fselc):
    # need to load psel first
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
    count = { 0: 0, 1: 0, 2: 0, }
    inv_type = flavor_inv_type
    idx = 0
    xg_list = psel.to_list()
    for xg, inv_acc in [ (xg, inv_acc) for xg in xg_list for inv_acc in (0, 1, 2,) ]:
        xg_str = f"({xg[0]},{xg[1]},{xg[2]},{xg[3]})"
        tag = f"xg={xg_str} ; type={inv_type} ; accuracy={inv_acc}"
        if not sfr.has(tag):
            continue
        q.displayln_info(0, f"load_prop_psrc_fsel: idx={idx} ; {tag} ; path_s={path_s}")
        idx += 1
        # load fsel psnk prop
        sc_prop = q.SelProp(fselc)
        sc_prop.load_double_from_float(sfr, tag)
        cache_fsel[f"{tag} ; psrc ; fsel"] = sc_prop
        # check psel psnk prop
        sp_prop_diff = q.PselProp(psel)
        sp_prop_diff @= cache_fsel[f"{tag} ; psrc ; fsel"]
        sp_prop_diff -= cache_psel[f"{tag} ; psrc ; psel"]
        assert sp_prop_diff.qnorm() <= 1e-14 * cache_psel[f"{tag} ; psrc ; psel"].qnorm()
        count[inv_acc] += 1
    sfr.close()
    assert cache_prob[f"type={flavor_inv_type} ; accuracy=0 ; psrc ; prob"] == count[0] / len(xg_list)
    assert cache_prob[f"type={flavor_inv_type} ; accuracy=1 ; psrc ; prob"] == rup.dict_params[job_tag]["prob_acc_1_psrc"]
    assert cache_prob[f"type={flavor_inv_type} ; accuracy=2 ; psrc ; prob"] == rup.dict_params[job_tag]["prob_acc_2_psrc"]

@q.timer
def load_prop_rand_u1_fsel(job_tag, traj, flavor, *, fsel):
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
    cache_fsel[f"type={inv_type} ; accuracy={inv_acc} ; rand_u1 ; fsel"] = s_prop_avg

### -------

@q.timer_verbose
def run_get_prop(job_tag, traj, *, get_gt, get_psel, get_fsel, get_psel_smear, get_wi):
    @q.timer_verbose
    def mk_get_prop():
        wi = get_wi()
        gt = get_gt()
        psel = get_psel()
        psel_smear = get_psel_smear()
        fsel, fselc = get_fsel()
        # ADJUST ME
        load_prop_wsrc_psel(job_tag, traj, "l", wi = wi, psel = psel, fsel = fsel, fselc = fselc, gt = gt)
        load_prop_wsrc_psel(job_tag, traj, "s", wi = wi, psel = psel, fsel = fsel, fselc = fselc, gt = gt)
        load_prop_wsrc_fsel(job_tag, traj, "l", wi = wi, psel = psel, fsel = fsel, fselc = fselc, gt = gt)
        load_prop_wsrc_fsel(job_tag, traj, "s", wi = wi, psel = psel, fsel = fsel, fselc = fselc, gt = gt)
        load_prop_psrc_psel(job_tag, traj, "l", psel = psel, fsel = fsel, fselc = fselc)
        load_prop_psrc_psel(job_tag, traj, "s", psel = psel, fsel = fsel, fselc = fselc)
        load_prop_psrc_fsel(job_tag, traj, "l", psel = psel, fsel = fsel, fselc = fselc)
        load_prop_psrc_fsel(job_tag, traj, "s", psel = psel, fsel = fsel, fselc = fselc)
        # load_prop_rand_u1_fsel(job_tag, traj, "l", fsel = fsel)
        # load_prop_rand_u1_fsel(job_tag, traj, "s", fsel = fsel)
        # load_prop_rand_u1_fsel(job_tag, traj, "c", fsel = fsel)
        #
        prop_cache = q.mk_cache(f"prop_cache", f"{job_tag}", f"{traj}")
        psel_pos_dict = dict([ (tuple(pos), i) for i, pos in enumerate(psel.to_list()) ])
        fsel_pos_dict = dict([ (tuple(pos), i) for i, pos in enumerate(fsel.to_psel_local().to_list()) ])
        fselc_pos_dict = dict([ (tuple(pos), i) for i, pos in enumerate(fselc.to_psel_local().to_list()) ])
        def get_prop(flavor, p_snk, p_src):
            return get_prop_snk_src(prop_cache, flavor, p_snk, p_src,
                    psel_pos_dict = psel_pos_dict,
                    fsel_pos_dict = fsel_pos_dict,
                    fselc_pos_dict = fselc_pos_dict)
        return get_prop
    return q.lazy_call(mk_get_prop)
