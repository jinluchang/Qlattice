from . import rbc_ukqcd as ru
from . import rbc_ukqcd_params as rup
from .jobs import *

# -----------------------------------------------------------------------------

@q.timer
def run_get_inverter(job_tag, traj, *, inv_type, get_gf, get_gt = None, get_eig = None):
    if None in [ get_gf, ]:
        return
    if get_gt is None:
        get_gt = lambda: None
    if get_eig is None:
        get_eig = lambda: None
    gf = get_gf()
    gt = get_gt()
    eig = get_eig()
    for inv_acc in [ 0, 1, 2, ]:
        ru.get_inv(gf, job_tag, inv_type, inv_acc, gt = gt, eig = eig)

# -----------------------------------------------------------------------------

@q.timer_verbose
def compute_prop_1(inv, src, *, tag, sfw, path_sp, psel, fsel, fselc):
    fn_sp = os.path.join(path_sp, f"{tag}.lat")
    fn_spw = os.path.join(path_sp, f"{tag} ; wsnk.lat")
    sol = inv * src
    s_sol = q.SelProp(fselc)
    s_sol @= sol
    s_sol.save_float_from_double(sfw, tag)
    sp_sol = q.PselProp(psel)
    sp_sol @= s_sol
    sp_sol.save(get_save_path(fn_sp))
    sfw.flush()
    sol_ws = sol.glb_sum_tslice()
    sol_ws.save(get_save_path(fn_spw))
    return sol

@q.timer
def compute_prop_wsrc(gf, gt, tslice, job_tag, inv_type, inv_acc, *,
        idx, sfw, path_sp, psel, fsel, fselc, eig, finished_tags):
    tag = f"tslice={tslice} ; type={inv_type} ; accuracy={inv_acc}"
    if tag in finished_tags:
        return None
    q.check_stop()
    q.check_time_limit()
    q.displayln_info(f"compute_prop_wsrc: idx={idx} tslice={tslice}", job_tag, inv_type, inv_acc)
    inv = ru.get_inv(gf, job_tag, inv_type, inv_acc, gt = gt, eig = eig)
    total_site = ru.get_total_site(job_tag)
    geo = q.Geometry(total_site, 1)
    src = q.mk_wall_src(geo, tslice)
    prop = compute_prop_1(inv, src, tag = tag, sfw = sfw, path_sp = path_sp,
                          psel = psel, fsel = fsel, fselc = fselc)

@q.timer_verbose
def compute_prop_wsrc_all(job_tag, traj, *,
                          inv_type, gf, gt, wi, psel, fsel, fselc, eig):
    inv_type_names = [ "light", "strange", ]
    inv_type_name = inv_type_names[inv_type]
    path_s = f"{job_tag}/prop-wsrc-{inv_type_name}/traj-{traj}"
    path_sp = f"{job_tag}/psel-prop-wsrc-{inv_type_name}/traj-{traj}"
    finished_tags = q.properly_truncate_fields(get_save_path(path_s + ".acc"))
    sfw = q.open_fields(get_save_path(path_s + ".acc"), "a", [ 1, 2, 2, 4, ])
    for inv_acc in [ 2, 1 ]:
        for p in wi:
            idx, tslice, inv_type_p, inv_acc_p = p
            if inv_type_p == inv_type and inv_acc_p == inv_acc:
                compute_prop_wsrc(gf, gt, tslice, job_tag, inv_type, inv_acc,
                        idx = idx, sfw = sfw, path_sp = path_sp,
                        psel = psel, fsel = fsel, fselc = fselc, eig = eig,
                        finished_tags = finished_tags)
    sfw.close()
    q.qtouch_info(get_save_path(os.path.join(path_sp, "checkpoint.txt")))
    # q.qtouch_info(get_save_path(os.path.join(path_sp, "checkpoint ; wsnk.txt")))
    q.qrename_info(get_save_path(path_s + ".acc"), get_save_path(path_s))
    q.qar_create_info(get_save_path(path_sp + ".qar"), get_save_path(path_sp), is_remove_folder_after = True)
    # q.qar_create_info(get_save_path(path_s + ".qar"), get_save_path(path_s), is_remove_folder_after = True)

@q.timer
def run_prop_wsrc(job_tag, traj, *, inv_type, get_gf, get_eig, get_gt, get_psel, get_fsel, get_wi):
    if None in [ get_gf, get_gt, get_psel, get_fsel, ]:
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
        fsel, fselc = get_fsel()
        psel = get_psel()
        wi = get_wi()
        compute_prop_wsrc_all(job_tag, traj,
                              inv_type = inv_type, gf = gf, gt = gt, wi = wi,
                              psel = psel, fsel = fsel, fselc = fselc, eig = eig)
        q.release_lock()

# -----------------------------------------------------------------------------

@q.timer_verbose
def compute_prop_2(inv, src, *, tag, sfw, path_sp, psel, fsel, fselc, gt):
    fn_sp = os.path.join(path_sp, f"{tag}.lat")
    fn_spw = os.path.join(path_sp, f"{tag} ; wsnk.lat")
    sol = inv * src
    s_sol = q.SelProp(fselc)
    s_sol @= sol
    s_sol.save_float_from_double(sfw, tag)
    sp_sol = q.PselProp(psel)
    sp_sol @= s_sol
    sp_sol.save(get_save_path(fn_sp))
    sfw.flush()
    sol_gt = gt * sol
    sol_ws = sol_gt.glb_sum_tslice()
    sol_ws.save(get_save_path(fn_spw))
    return sol

@q.timer
def compute_prop_psrc(job_tag, xg_src, inv_type, inv_acc, *,
        idx, gf, gt, sfw, path_sp, psel, fsel, fselc, eig, finished_tags):
    xg = xg_src
    xg_str = f"({xg[0]},{xg[1]},{xg[2]},{xg[3]})"
    tag = f"xg={xg_str} ; type={inv_type} ; accuracy={inv_acc}"
    if tag in finished_tags:
        return None
    q.check_stop()
    q.check_time_limit()
    q.displayln_info(f"compute_prop_psrc: {job_tag} idx={idx} tag='{tag}'")
    inv = ru.get_inv(gf, job_tag, inv_type, inv_acc, eig = eig)
    total_site = get_param(job_tag, "total_site")
    geo = q.Geometry(total_site, 1)
    src = q.mk_point_src(geo, xg_src)
    prop = compute_prop_2(inv, src, tag = tag, sfw = sfw, path_sp = path_sp,
                          psel = psel, fsel = fsel, fselc = fselc, gt = gt)

@q.timer_verbose
def compute_prop_psrc_all(job_tag, traj, *,
                          inv_type, gf, gt, psel, fsel, fselc, eig):
    inv_type_names = [ "light", "strange", ]
    inv_type_name = inv_type_names[inv_type]
    path_s = f"{job_tag}/prop-psrc-{inv_type_name}/traj-{traj}"
    path_sp = f"{job_tag}/psel-prop-psrc-{inv_type_name}/traj-{traj}"
    finished_tags = q.properly_truncate_fields(get_save_path(path_s + ".acc"))
    sfw = q.open_fields(get_save_path(path_s + ".acc"), "a", [ 1, 2, 2, 4, ])
    def comp(idx, xg_src, inv_acc):
        compute_prop_psrc(job_tag, xg_src, inv_type, inv_acc,
                idx = idx, gf = gf, gt = gt, sfw = sfw, path_sp = path_sp,
                psel = psel, fsel = fsel, fselc = fselc,
                eig = eig, finished_tags = finished_tags)
    prob1 = get_param(job_tag, "prob_acc_1_psrc")
    prob2 = get_param(job_tag, "prob_acc_2_psrc")
    rs = q.RngState(f"seed {job_tag} {traj}").split(f"compute_prop_psrc_all(ama)")
    for idx, xg_src in enumerate(psel.to_list()):
        r = rs.split(f"{tuple(xg_src)}").u_rand_gen()
        assert 0 <= r and r <= 1
        comp(idx, xg_src, inv_acc = 0)
        if r <= prob1:
            comp(idx, xg_src, inv_acc = 1)
        if r <= prob2:
            comp(idx, xg_src, inv_acc = 2)
    sfw.close()
    q.qtouch_info(get_save_path(os.path.join(path_sp, "checkpoint.txt")))
    q.qrename_info(get_save_path(path_s + ".acc"), get_save_path(path_s))
    q.qar_create_info(get_save_path(path_sp + ".qar"), get_save_path(path_sp), is_remove_folder_after = True)
    # q.qar_create_info(get_save_path(path_s + ".qar"), get_save_path(path_s), is_remove_folder_after = True)

@q.timer
def run_prop_psrc(job_tag, traj, *, inv_type, get_gf, get_eig, get_gt, get_psel, get_fsel):
    if None in [ get_gf, get_gt, get_psel, get_fsel, ]:
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
        fsel, fselc = get_fsel()
        psel = get_psel()
        compute_prop_psrc_all(job_tag, traj,
                              inv_type = inv_type, gf = gf, gt = gt,
                              psel = psel, fsel = fsel, fselc = fselc, eig = eig)
        q.release_lock()

# -----------------------------------------------------------------------------

@q.timer_verbose
def compute_prop_rand_u1_type_acc(*, sfw, job_tag, traj, gf, eig, fsel, idx_rand_u1, inv_type, inv_acc, finished_tags):
    # same rand source for different inv_type
    tag = f"idx_rand_u1={idx_rand_u1} ; type={inv_type} ; accuracy={inv_acc}"
    if tag in finished_tags:
        return
    q.check_stop()
    q.check_time_limit()
    inv = ru.get_inv(gf, job_tag, inv_type, inv_acc, eig = eig)
    rs = q.RngState(f"seed {job_tag} {traj}").split(f"compute_prop_rand_u1(rand_u1)").split(str(idx_rand_u1))
    s_prop = q.mk_rand_u1_prop(inv, fsel, rs)
    s_prop.save_float_from_double(sfw, tag)
    sfw.flush()
    return s_prop

@q.timer_verbose
def compute_prop_rand_u1(*, job_tag, traj, inv_type, gf, path_s, fsel, eig = None):
    # use fsel instead of fselc
    n_rand_u1_fsel = get_param(job_tag, "n_rand_u1_fsel")
    total_site = get_param(job_tag, "total_site")
    geo = q.Geometry(total_site, 1)
    finished_tags = q.properly_truncate_fields(get_save_path(path_s + ".acc"))
    sfw = q.open_fields(get_save_path(path_s + ".acc"), "a", [ 1, 2, 2, 4, ])
    def comp(idx_rand_u1, inv_acc):
        compute_prop_rand_u1_type_acc(
                sfw = sfw,
                job_tag = job_tag, traj = traj,
                gf = gf, eig = eig, fsel = fsel,
                idx_rand_u1 = idx_rand_u1,
                inv_type = inv_type, inv_acc = inv_acc,
                finished_tags = finished_tags,)
    rs = q.RngState(f"seed {job_tag} {traj}").split(f"compute_prop_rand_u1(ama)")
    prob1 = get_param(job_tag, "prob_acc_1_rand_u1")
    prob2 = get_param(job_tag, "prob_acc_2_rand_u1")
    for idx_rand_u1 in range(n_rand_u1_fsel):
        r = rs.split(str(idx_rand_u1)).u_rand_gen()
        inv_acc = 0
        assert 0 <= r and r <= 1
        comp(idx_rand_u1, inv_acc)
        inv_acc = 1
        if r <= prob1:
            comp(idx_rand_u1, inv_acc)
        inv_acc = 2
        if r <= prob2:
            comp(idx_rand_u1, inv_acc)
    sfw.close()
    q.qrename_info(get_save_path(path_s + ".acc"), get_save_path(path_s))
    # q.qar_create_info(get_save_path(path_s + ".qar"), get_save_path(path_s), is_remove_folder_after = True)

@q.timer_verbose
def run_prop_rand_u1(job_tag, traj, *, inv_type, get_gf, get_fsel, get_eig = None):
    if None in [ get_gf, get_fsel, ]:
        return
    if get_eig is None:
        if inv_type == 0:
            return
        get_eig = lambda: None
    inv_type_names = [ "light", "strange", "charm", ]
    inv_type_name = inv_type_names[inv_type]
    path_s = f"{job_tag}/prop-rand-u1-{inv_type_name}/traj-{traj}"
    if get_load_path(path_s + "/geon-info.txt") is not None:
        return
    if q.obtain_lock(f"locks/{job_tag}-{traj}-rand-u1-{inv_type_name}"):
        gf = get_gf()
        fsel, fselc = get_fsel()
        eig = get_eig()
        compute_prop_rand_u1(
                job_tag = job_tag, traj = traj,
                inv_type = inv_type,
                gf = gf,
                path_s = path_s,
                fsel = fsel,
                eig = eig)
        q.release_lock()

# -----------------------------------------------------------------------------

@q.timer_verbose
def compute_prop_3(inv, src_smear, *, tag, sfw, path_sp, psel, fsel, fselc, gt, psel_smear, smear):
    fn_sp = os.path.join(path_sp, f"{tag}.lat")
    fn_spw = os.path.join(path_sp, f"{tag} ; wsnk.lat")
    fn_sps = os.path.join(path_sp, f"{tag} ; smear-snk.lat")
    sol = inv * src_smear
    s_sol = q.SelProp(fselc)
    s_sol @= sol
    s_sol.save_float_from_double(sfw, tag)
    sp_sol = q.PselProp(psel)
    sp_sol @= s_sol
    sp_sol.save(get_save_path(fn_sp))
    sfw.flush()
    sol_gt = gt * sol
    sol_ws = sol_gt.glb_sum_tslice()
    sol_ws.save(get_save_path(fn_spw))
    sol_smear = sol.copy()
    smear(sol_smear)
    sol_smear_psel = q.PselProp(psel_smear)
    sol_smear_psel @= sol_smear
    sol_smear_psel.save(get_save_path(fn_sps))
    return sol

@q.timer
def compute_prop_smear(job_tag, xg_src, inv_type, inv_acc, *,
        idx, gf, gt, sfw, path_sp, psel, fsel, fselc, psel_smear, gf_ape, eig, finished_tags):
    xg = xg_src
    xg_str = f"({xg[0]},{xg[1]},{xg[2]},{xg[3]})"
    tag = f"smear ; xg={xg_str} ; type={inv_type} ; accuracy={inv_acc}"
    if tag in finished_tags:
        return None
    q.check_stop()
    q.check_time_limit()
    q.displayln_info(f"compute_prop_smear: {job_tag} idx={idx} tag='{tag}'")
    inv = ru.get_inv(gf, job_tag, inv_type, inv_acc, eig = eig)
    total_site = ru.get_total_site(job_tag)
    geo = q.Geometry(total_site, 1)
    coef = get_param(job_tag, "prop_smear_coef")
    step = get_param(job_tag, "prop_smear_step")
    def smear(src):
        q.prop_smear(src, gf_ape, coef, step)
    src = q.mk_point_src(geo, xg_src)
    smear(src)
    prop = compute_prop_3(inv, src, tag = tag, sfw = sfw, path_sp = path_sp,
                          psel = psel, fsel = fsel, fselc = fselc, gt = gt,
                          psel_smear = psel_smear, smear = smear)

@q.timer_verbose
def compute_prop_smear_all(job_tag, traj, *,
        inv_type, gf, gt, psel, fsel, fselc, psel_smear, gf_ape, eig,
        ):
    inv_type_names = [ "light", "strange", ]
    inv_type_name = inv_type_names[inv_type]
    path_s = f"{job_tag}/prop-smear-{inv_type_name}/traj-{traj}"
    path_sp = f"{job_tag}/psel-prop-smear-{inv_type_name}/traj-{traj}"
    finished_tags = q.properly_truncate_fields(get_save_path(path_s + ".acc"))
    sfw = q.open_fields(get_save_path(path_s + ".acc"), "a", [ 1, 2, 2, 4, ])
    def comp(idx, xg_src, inv_acc):
        compute_prop_smear(job_tag, xg_src, inv_type, inv_acc,
                idx = idx, gf = gf, gt = gt, sfw = sfw, path_sp = path_sp,
                psel = psel, fsel = fsel, fselc = fselc,
                psel_smear = psel_smear, gf_ape = gf_ape,
                eig = eig, finished_tags = finished_tags)
    prob1 = get_param(job_tag, "prob_acc_1_smear")
    prob2 = get_param(job_tag, "prob_acc_2_smear")
    rs = q.RngState(f"seed {job_tag} {traj}").split(f"compute_prop_smear_all(ama)")
    for idx, xg_src in enumerate(psel_smear.to_list()):
        r = rs.split(f"{tuple(xg_src)}").u_rand_gen()
        assert 0 <= r and r <= 1
        comp(idx, xg_src, inv_acc = 0)
        if r <= prob1:
            comp(idx, xg_src, inv_acc = 1)
        if r <= prob2:
            comp(idx, xg_src, inv_acc = 2)
    sfw.close()
    q.qtouch_info(get_save_path(os.path.join(path_sp, "checkpoint.txt")))
    q.qrename_info(get_save_path(path_s + ".acc"), get_save_path(path_s))
    q.qar_create_info(get_save_path(path_sp + ".qar"), get_save_path(path_sp), is_remove_folder_after = True)
    # q.qar_create_info(get_save_path(path_s + ".qar"), get_save_path(path_s), is_remove_folder_after = True)

@q.timer
def run_prop_smear(job_tag, traj, *, inv_type, get_gf, get_gf_ape, get_eig, get_gt, get_psel, get_fsel, get_psel_smear):
    if None in [ get_gf, get_gt, get_gf_ape, get_psel, get_fsel, ]:
        return
    if get_eig is None:
        if inv_type == 0:
            return
        get_eig = lambda: None
    inv_type_names = [ "light", "strange", ]
    inv_type_name = inv_type_names[inv_type]
    if get_load_path(f"{job_tag}/prop-smear-{inv_type_name}/traj-{traj}/geon-info.txt") is not None:
        return
    if q.obtain_lock(f"locks/{job_tag}-{traj}-smear-{inv_type_name}"):
        gf = get_gf()
        gt = get_gt()
        eig = get_eig()
        psel = get_psel()
        fsel, fselc = get_fsel()
        psel_smear = get_psel_smear()
        gf_ape = get_gf_ape()
        compute_prop_smear_all(job_tag, traj,
                inv_type = inv_type, gf = gf, gf_ape = gf_ape, gt = gt,
                psel = psel, fsel = fsel, fselc = fselc, eig = eig, psel_smear = psel_smear)
        q.release_lock()

# -----------------------------------------------------------------------------
