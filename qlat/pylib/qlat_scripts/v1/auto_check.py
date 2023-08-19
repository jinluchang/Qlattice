from .jobs import *

# ----

@q.timer
def get_all_points(total_site):
    n_points = total_site[0] * total_site[1] * total_site[2] * total_site[3]
    xg_list = []
    for index in range(n_points):
        xg = q.Coordinate()
        xg.from_index(index, total_site)
        xg_list.append(xg)
    return xg_list

@q.timer
def get_all_points_psel(total_site):
    geo = q.Geometry(total_site, 1)
    xg_list = get_all_points(total_site)
    psel = q.PointsSelection([ xg.to_list() for xg in xg_list ], geo)
    return psel

# ----

@q.timer
def run_get_inverter_checker(job_tag, traj, *, inv_type, get_gf, get_gt = None, get_eig = None):
    if None in [ get_gf, ]:
        return
    if get_gt is None:
        get_gt = lambda: None
    if get_eig is None:
        get_eig = lambda: None
    gf = get_gf()
    gt = get_gt()
    eig = get_eig()
    inv_acc = 2
    from . import rbc_ukqcd as ru
    ru.get_inv(gf, job_tag, inv_type, inv_acc, gt = gt, eig = eig)

@q.timer_verbose
def compute_prop_1_checker(inv, src, *, tag, sfw, path_sp):
    fn_spw = os.path.join(path_sp, f"{tag} ; wsnk.lat")
    sol = inv * src
    sol.save_float_from_double(sfw, tag)
    sfw.flush()
    sol_ws = sol.glb_sum_tslice()
    sol_ws.save(get_save_path(fn_spw))
    return sol

@q.timer
def compute_prop_wsrc_checker(job_tag, tslice, inv_type, inv_acc, *,
                              idx, gf, gt, sfw, path_sp, eig, finished_tags):
    from . import rbc_ukqcd as ru
    tag = f"tslice={tslice} ; type={inv_type} ; accuracy={inv_acc}"
    if tag in finished_tags:
        return None
    q.displayln_info(f"compute_prop_wsrc: idx={idx} tslice={tslice}", job_tag, inv_type, inv_acc)
    inv = ru.get_inv(gf, job_tag, inv_type, inv_acc, gt = gt, eig = eig)
    total_site = q.Coordinate(get_param(job_tag, "total_site"))
    geo = q.Geometry(total_site, 1)
    src = q.mk_wall_src(geo, tslice)
    prop = compute_prop_1_checker(inv, src, tag = tag, sfw = sfw, path_sp = path_sp)

@q.timer_verbose
def compute_prop_wsrc_all_checker(job_tag, traj, *,
                                  inv_type, gf, gt, eig):
    total_site = q.Coordinate(get_param(job_tag, "total_site"))
    inv_type_names = [ "light", "strange", ]
    inv_type_name = inv_type_names[inv_type]
    path_s = f"{job_tag}/prop-wsrc-{inv_type_name}/traj-{traj}"
    path_sp = f"{job_tag}/psel-prop-wsrc-{inv_type_name}/traj-{traj}"
    finished_tags = q.properly_truncate_fields(get_save_path(path_s + ".acc"))
    sfw = q.open_fields(get_save_path(path_s + ".acc"), "a", [ 1, 1, 1, 4, ])
    inv_acc = 2
    for idx, tslice in enumerate(range(total_site[3])):
        compute_prop_wsrc_checker(job_tag, tslice, inv_type, inv_acc = 2,
                                  idx = idx, gf = gf, gt = gt, sfw = sfw, path_sp = path_sp,
                                  eig = eig, finished_tags = finished_tags)
    sfw.close()
    q.qtouch_info(get_save_path(os.path.join(path_sp, "checkpoint.txt")))
    q.qrename_info(get_save_path(path_s + ".acc"), get_save_path(path_s))
    q.qar_create_info(get_save_path(path_sp + ".qar"), get_save_path(path_sp), is_remove_folder_after = True)
    q.qar_create_info(get_save_path(path_s + ".qar"), get_save_path(path_s), is_remove_folder_after = True)

@q.timer
def run_prop_wsrc_checker(job_tag, traj, *, inv_type, get_gf, get_eig, get_gt):
    if None in [ get_gf, get_gt, ]:
        return
    if get_eig is None:
        if inv_type == 0:
            return
        get_eig = lambda: None
    inv_type_names = [ "light", "strange", ]
    inv_type_name = inv_type_names[inv_type]
    if get_load_path(f"{job_tag}/prop-wsrc-{inv_type_name}/traj-{traj}/geon-info.txt") is not None:
        return
    if q.obtain_lock(f"locks/{job_tag}-{traj}-wsrc-{inv_type_name}"):
        gf = get_gf()
        gt = get_gt()
        eig = get_eig()
        compute_prop_wsrc_all_checker(job_tag, traj,
                                      inv_type = inv_type, gf = gf, gt = gt, eig = eig)
        q.release_lock()

# ----

@q.timer_verbose
def compute_prop_2_checker(inv, src, *, tag, sfw):
    sol = inv * src
    sol.save_float_from_double(sfw, tag)
    sfw.flush()
    return sol

@q.timer
def compute_prop_psrc_checker(job_tag, xg_src, inv_type, inv_acc, *,
                              idx, gf, gt, sfw, eig, finished_tags):
    from . import rbc_ukqcd as ru
    xg = xg_src.to_list()
    xg_str = f"({xg[0]},{xg[1]},{xg[2]},{xg[3]})"
    tag = f"xg={xg_str} ; type={inv_type} ; accuracy={inv_acc}"
    if tag in finished_tags:
        return None
    q.displayln_info(f"compute_prop_psrc: {job_tag} idx={idx} tag='{tag}'")
    inv = ru.get_inv(gf, job_tag, inv_type, inv_acc, gt = gt, eig = eig)
    total_site = q.Coordinate(get_param(job_tag, "total_site"))
    geo = q.Geometry(total_site, 1)
    src = q.mk_point_src(geo, xg)
    prop = compute_prop_2_checker(inv, src, tag = tag, sfw = sfw)

@q.timer_verbose
def compute_prop_psrc_all_checker(job_tag, traj, *,
                                  inv_type, gf, gt, eig):
    total_site = q.Coordinate(get_param(job_tag, "total_site"))
    inv_type_names = [ "light", "strange", ]
    inv_type_name = inv_type_names[inv_type]
    path_s = f"{job_tag}/prop-psrc-{inv_type_name}/traj-{traj}"
    finished_tags = q.properly_truncate_fields(get_save_path(path_s + ".acc"))
    sfw = q.open_fields(get_save_path(path_s + ".acc"), "a", [ 1, 1, 1, 4, ])
    def comp(idx, xg_src, inv_acc):
        compute_prop_psrc_checker(job_tag, xg_src, inv_type, inv_acc,
                idx = idx, gf = gf, gt = gt, sfw = sfw,
                eig = eig, finished_tags = finished_tags)
    for idx, xg_src in enumerate(get_all_points(total_site)):
        comp(idx, xg_src, inv_acc = 2)
    sfw.close()
    q.qrename_info(get_save_path(path_s + ".acc"), get_save_path(path_s))
    q.qar_create_info(get_save_path(path_s + ".qar"), get_save_path(path_s), is_remove_folder_after = True)

@q.timer
def run_prop_psrc_checker(job_tag, traj, *, inv_type, get_gf, get_eig, get_gt):
    if None in [ get_gf, get_gt, ]:
        return
    if get_eig is None:
        if inv_type == 0:
            return
        get_eig = lambda: None
    inv_type_names = [ "light", "strange", ]
    inv_type_name = inv_type_names[inv_type]
    if get_load_path(f"{job_tag}/prop-psrc-{inv_type_name}/traj-{traj}/geon-info.txt") is not None:
        return
    if q.obtain_lock(f"locks/{job_tag}-{traj}-psrc-{inv_type_name}"):
        gf = get_gf()
        gt = get_gt()
        eig = get_eig()
        compute_prop_psrc_all_checker(job_tag, traj,
                                      inv_type = inv_type, gf = gf, gt = gt,
                                      eig = eig)
        q.release_lock()

# ----

@q.timer_verbose
def load_prop_psrc(job_tag, traj, inv_type):
    inv_tag_list = [ "l", "s", ]
    inv_tag = inv_tag_list[inv_type]
    cache = q.mk_cache(f"prop_cache", f"{job_tag}", f"{traj}", inv_tag)
    if "psnk-psrc" in cache:
        return
    total_site = q.Coordinate(get_param(job_tag, "total_site"))
    inv_type_names = [ "light", "strange", ]
    inv_type_name = inv_type_names[inv_type]
    inv_acc = 2
    path_s = f"{job_tag}/prop-psrc-{inv_type_name}/traj-{traj}/geon-info.txt"
    psel = get_all_points_psel(total_site)
    prop_list = []
    xg_list = [ q.Coordinate(xg) for xg in psel.to_list() ]
    sfr = q.open_fields(get_load_path(path_s), "r")
    for xg_src in xg_list:
        xg_idx = xg_src.to_index(total_site)
        xg = xg_src.to_list()
        xg_str = f"({xg[0]},{xg[1]},{xg[2]},{xg[3]})"
        tag = f"xg={xg_str} ; type={inv_type} ; accuracy={inv_acc}"
        prop = q.Prop()
        prop.load_double_from_float(sfr, tag)
        sp_prop = q.PselProp(psel)
        sp_prop @= prop
        prop_list.append(sp_prop)
    sfr.close()
    cache["psnk-psrc"] = prop_list

@q.timer_verbose
def load_prop_wsrc(job_tag, traj, inv_type):
    inv_tag_list = [ "l", "s", ]
    inv_tag = inv_tag_list[inv_type]
    cache = q.mk_cache(f"prop_cache", f"{job_tag}", f"{traj}", inv_tag)
    if "psnk-wsrc" in cache and "wsnk-wsrc" in cache:
        return
    total_site = q.Coordinate(get_param(job_tag, "total_site"))
    inv_type_names = [ "light", "strange", ]
    inv_type_name = inv_type_names[inv_type]
    inv_acc = 2
    path_s = f"{job_tag}/prop-wsrc-{inv_type_name}/traj-{traj}/geon-info.txt"
    path_sp = f"{job_tag}/psel-prop-wsrc-{inv_type_name}/traj-{traj}/"
    psel = get_all_points_psel(total_site)
    psel_ts = q.get_psel_tslice(total_site)
    prop_list = []
    prop2_list = []
    tslice_list = list(range(total_site[3]))
    sfr = q.open_fields(get_load_path(path_s), "r")
    for tslice in tslice_list:
        tag = f"tslice={tslice} ; type={inv_type} ; accuracy={inv_acc}"
        prop = q.Prop()
        prop.load_double_from_float(sfr, tag)
        sp_prop = q.PselProp(psel)
        sp_prop @= prop
        prop_list.append(sp_prop)
        fn_spw = os.path.join(path_sp, f"{tag} ; wsnk.lat")
        sp_prop2 = q.PselProp(psel_ts)
        sp_prop2.load(get_load_path(fn_spw))
        prop2_list.append(sp_prop2)
    sfr.close()
    cache["psnk-wsrc"] = prop_list
    cache["wsnk-wsrc"] = prop2_list

@q.timer_verbose
def run_get_prop_checker(job_tag, traj, *,
                         get_gf,
                         get_gt):
    traj_gf = traj
    fns_props = [
            (f"{job_tag}/prop-psrc-light/traj-{traj_gf}.qar", f"{job_tag}/prop-psrc-light/traj-{traj_gf}/geon-info.txt",),
            (f"{job_tag}/prop-psrc-strange/traj-{traj_gf}.qar", f"{job_tag}/prop-psrc-strange/traj-{traj_gf}/geon-info.txt",),
            (f"{job_tag}/prop-wsrc-light/traj-{traj_gf}.qar", f"{job_tag}/prop-wsrc-light/traj-{traj_gf}/geon-info.txt",),
            (f"{job_tag}/prop-wsrc-strange/traj-{traj_gf}.qar", f"{job_tag}/prop-wsrc-strange/traj-{traj_gf}/geon-info.txt",),
            (f"{job_tag}/psel-prop-wsrc-light/traj-{traj_gf}.qar", f"{job_tag}/psel-prop-wsrc-light/traj-{traj_gf}/checkpoint.txt",),
            (f"{job_tag}/psel-prop-wsrc-strange/traj-{traj_gf}.qar", f"{job_tag}/psel-prop-wsrc-strange/traj-{traj_gf}/checkpoint.txt",),
            ]
    for fn in fns_props:
        if get_load_path(fn) is None:
            return None
    @q.timer_verbose
    def mk_get_prop():
        q.timer_fork()
        total_site = q.Coordinate(get_param(job_tag, "total_site"))
        gf = get_gf()
        gt = get_gt()
        #
        load_prop_psrc(job_tag, traj, inv_type = 0)
        load_prop_psrc(job_tag, traj, inv_type = 1)
        load_prop_wsrc(job_tag, traj, inv_type = 0)
        load_prop_wsrc(job_tag, traj, inv_type = 1)
        #
        prop_cache = q.mk_cache(f"prop_cache", f"{job_tag}", f"{traj}")
        def get_prop(flavor, p_snk, p_src):
            cache = prop_cache[flavor]
            p_snk_tag, p_snk_xg = p_snk
            p_src_tag, p_src_xg = p_src
            if p_snk_tag == "point" and p_src_tag == "point":
                prop_list = cache["psnk-psrc"]
                p_src_idx = p_src_xg.to_index(total_site)
                p_snk_idx = p_snk_xg.to_index(total_site)
                return prop_list[p_src_idx].get_elem_wm(p_snk_idx)
            elif p_snk_tag == "point" and p_src_tag == "wall":
                prop_list = cache["psnk-wsrc"]
                p_snk_idx = p_snk_xg.to_index(total_site)
                return prop_list[p_src_xg].get_elem_wm(p_snk_idx)
            elif p_snk_tag == "wall" and p_src_tag == "point":
                prop_list = cache["psnk-wsrc"]
                p_src_idx = p_src_xg.to_index(total_site)
                return q.wilson_matrix_g5_herm(prop_list[p_snk_xg].get_elem_wm(p_src_idx))
            elif p_snk_tag == "wall" and p_src_tag == "wall":
                prop_list = cache["wsnk-wsrc"]
                return prop_list[p_src_xg].get_elem_wm(p_snk_xg)
            else:
                raise Exception(f"get_prop: f={flavor} snk={p_snk} src={p_src}")
        #
        q.timer_display()
        q.timer_merge()
        return get_prop
    return q.lazy_call(mk_get_prop)

# ----
