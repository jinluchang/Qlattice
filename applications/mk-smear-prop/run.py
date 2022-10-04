#!/usr/bin/env python3

# Need --mpi X.X.X.X --mpi X.X.X runtime option

import qlat_gpt as qg
import rbc_ukqcd as ru

from jobs import *

import params

load_path_list[:] = [
        "results",
        "../mk-gf-gt/results",
        "../mk-sel/results",
        "../mk-lanc/results",
        "../qcddata",
        "../qcddata-1",
        "../qcddata-2",
        "../qcddata-3",
        "../qcddata-4",
        os.path.join(os.getenv("HOME"), "Qlat-sample-data/mk-gf-gt/results"),
        os.path.join(os.getenv("HOME"), "Qlat-sample-data/mk-sel/results"),
        os.path.join(os.getenv("HOME"), "Qlat-sample-data/mk-lanc/results"),
        os.path.join(os.getenv("HOME"), "qcddata"),
        ]

@q.timer_verbose
def compute_prop(inv, xg_src, *, job_tag, sfw, tag, path_sp, psel, fsel, fselc, psel_smear, gf_ape, gt):
    fn_sp = os.path.join(path_sp, f"{tag}.lat")
    fn_spw = os.path.join(path_sp, f"{tag} ; wsnk.lat")
    fn_sps = os.path.join(path_sp, f"{tag} ; smear-snk.lat")
    total_site = ru.get_total_site(job_tag)
    geo = q.Geometry(total_site, 1)
    src = q.mk_point_src(geo, xg_src)
    coef = rup.dict_params[job_tag]["prop_smear_coef"]
    step = rup.dict_params[job_tag]["prop_smear_step"]
    src_smear = src.copy()
    q.prop_smear(src_smear, gf_ape, coef, step)
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
    q.prop_smear(sol_smear, gf_ape, coef, step)
    sol_smear_psel = q.PselProp(psel_smear)
    sol_smear_psel @= sol_smear
    sol_smear_psel.save(get_save_path(fn_sps))

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
    compute_prop(inv, xg_src, job_tag = job_tag, sfw = sfw, tag = tag, path_sp = path_sp,
            psel = psel, fsel = fsel, fselc = fselc,
            psel_smear = psel_smear, gf_ape = gf_ape, gt = gt)

@q.timer_verbose
def compute_prop_smear_all(job_tag, traj, *,
        inv_type, gf, gt, psel, fsel, fselc, psel_smear, gf_ape, eig,
        ):
    inv_type_names = [ "light", "strange", ]
    inv_type_name = inv_type_names[inv_type]
    path_s = f"{job_tag}/prop-smear-{inv_type_name}/traj-{traj}"
    path_sp = f"{job_tag}/psel-prop-smear-{inv_type_name}/traj-{traj}"
    finished_tags = q.properly_truncate_fields(get_save_path(path_s + ".acc"))
    sfw = q.open_fields(get_save_path(path_s + ".acc"), "a", [ 1, 1, 1, 4, ])
    def comp(idx, xg_src, inv_acc):
        compute_prop_smear(job_tag, xg_src, inv_type, inv_acc,
                idx = idx, gf = gf, gt = gt, sfw = sfw, path_sp = path_sp,
                psel = psel, fsel = fsel, fselc = fselc,
                psel_smear = psel_smear, gf_ape = gf_ape,
                eig = eig, finished_tags = finished_tags)
    prob1 = rup.dict_params[job_tag]["prob_acc_1_smear"]
    prob2 = rup.dict_params[job_tag]["prob_acc_2_smear"]
    rs = q.RngState(f"seed {job_tag} {traj}").split(f"compute_prop_smear_all(ama)")
    for idx, xg_src in enumerate(psel_smear.to_list()):
        r = rs.split(f"{tuple(xg_src)}").u_rand_gen()
        assert 0 <= r and r <= 1
        comp(idx, xg_src, inv_acc = 0)
        if r <= prob1:
            comp(idx, xg_src, inv_acc = 1)
        if r <= prob2:
            comp(idx, xg_src, inv_acc = 2)
    q.clean_cache(q.cache_inv)
    sfw.close()
    q.qtouch_info(get_save_path(os.path.join(path_sp, "checkpoint.txt")))
    q.qrename_info(get_save_path(path_s + ".acc"), get_save_path(path_s))
    q.qar_create_info(get_save_path(path_sp + ".qar"), get_save_path(path_sp), is_remove_folder_after = True)
    q.qar_create_info(get_save_path(path_s + ".qar"), get_save_path(path_s), is_remove_folder_after = True)

@q.timer
def run_prop_smear(job_tag, traj, *, inv_type, get_gf, get_gf_ape, get_eig, get_gt, get_psel, get_fsel, get_psel_smear):
    if None in [ get_gf, get_gt, get_gf_ape, get_eig, get_psel, get_fsel, ]:
        return
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

@q.timer_verbose
def run_job(job_tag, traj):
    fns_produce = [
            f"{job_tag}/prop-smear-light/traj-{traj}/geon-info.txt",
            f"{job_tag}/psel-prop-smear-light/traj-{traj}/checkpoint.txt",
            f"{job_tag}/prop-smear-strange/traj-{traj}/geon-info.txt",
            f"{job_tag}/psel-prop-smear-strange/traj-{traj}/checkpoint.txt",
            ]
    fns_need = [
            (f"{job_tag}/configs/ckpoint_lat.{traj}", f"{job_tag}/configs/ckpoint_lat.IEEE64BIG.{traj}",),
            f"{job_tag}/gauge-transform/traj-{traj}.field",
            f"{job_tag}/point-selection/traj-{traj}.txt",
            f"{job_tag}/field-selection/traj-{traj}.field",
            f"{job_tag}/point-selection-smear/traj-{traj}.txt",
            f"{job_tag}/eig/traj-{traj}",
            f"{job_tag}/eig/traj-{traj}/metadata.txt",
            f"{job_tag}/eig/traj-{traj}/eigen-values.txt",
            f"{job_tag}/eig-strange/traj-{traj}",
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
    assert get_psel is not None
    assert get_fsel is not None
    #
    get_gf_ape = run_gf_ape(job_tag, get_gf)
    #
    get_psel_smear = run_psel_smear(job_tag, traj)
    assert get_psel_smear is not None
    #
    def run_prop(inv_type, get_eig):
        run_prop_smear(job_tag, traj,
                inv_type = inv_type,
                get_gf = get_gf,
                get_gf_ape = get_gf_ape,
                get_eig = get_eig,
                get_gt = get_gt,
                get_psel = get_psel,
                get_fsel = get_fsel,
                get_psel_smear = get_psel_smear,
                )
    #
    def run_with_eig():
        get_eig = run_eig(job_tag, traj_gf, get_gf)
        run_prop(inv_type = 0, get_eig = get_eig)
    #
    def run_with_eig_strange():
        get_eig_strange = run_eig_strange(job_tag, traj_gf, get_gf)
        # get_eig_strange = lambda : None
        run_prop(inv_type = 1, get_eig = get_eig_strange)
    #
    run_with_eig()
    run_with_eig_strange()
    #
    q.clean_cache()
    q.timer_display()

qg.begin_with_gpt()

# ADJUST ME
job_tags = [
        "test-4nt8", "test-4nt16",
        # "32IH1",
        # "32IH2",
        # "24IH1",
        # "24IH2",
        # "24IH3",
        # "24D",
        # "32Dfine",
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

qg.end_with_gpt()
