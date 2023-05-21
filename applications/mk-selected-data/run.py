#!/usr/bin/env python3

# Need --mpi X.X.X.X --mpi X.X.X runtime option

import qlat_gpt as qg
import qlat_scripts.v1.rbc_ukqcd as ru

from qlat_scripts.v1.jobs import *
from qlat_scripts.v1.gen_data import *

load_path_list[:] = [
        "results",
        "qcddata",
        "qcddata-1",
        "qcddata-2",
        "qcddata-3",
        "qcddata-4",
        "qcddata-5",
        "../qcddata",
        "../qcddata-1",
        "../qcddata-2",
        "../qcddata-3",
        "../qcddata-4",
        "../qcddata-5",
        "../mk-gf-gt/results",
        "../mk-sel/results",
        "../mk-lanc/results",
        os.path.join(os.getenv("HOME"), "qcddata"),
        os.path.join(os.getenv("HOME"), "qcddata-1"),
        os.path.join(os.getenv("HOME"), "qcddata-2"),
        os.path.join(os.getenv("HOME"), "qcddata-3"),
        os.path.join(os.getenv("HOME"), "qcddata-4"),
        os.path.join(os.getenv("HOME"), "qcddata-5"),
        os.path.join(os.getenv("HOME"), "Qlat-sample-data/default/mk-gf-gt/results"),
        os.path.join(os.getenv("HOME"), "Qlat-sample-data/default/mk-sel/results"),
        os.path.join(os.getenv("HOME"), "Qlat-sample-data/default/mk-lanc/results"),
        ]

@q.timer_verbose
def compute_prop(inv, src, *, tag, sfw, fn_sp, psel, fsel, fselc):
    sol = inv * src
    s_sol = q.SelProp(fselc)
    s_sol @= sol
    s_sol.save_float_from_double(sfw, tag)
    sp_sol = q.PselProp(psel)
    sp_sol @= s_sol
    sp_sol.save(get_save_path(fn_sp))
    sfw.flush()
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
    fn_sp = os.path.join(path_sp, f"{tag}.lat")
    prop = compute_prop(inv, src, tag = tag, sfw = sfw, fn_sp = fn_sp, psel = psel, fsel = fsel, fselc = fselc)
    fn_spw = os.path.join(path_sp, f"{tag} ; wsnk.lat")
    prop.glb_sum_tslice().save(get_save_path(fn_spw))

@q.timer_verbose
def compute_prop_wsrc_all(gf, gt, wi, job_tag, inv_type, *,
        path_s, path_sp, psel, fsel, fselc, eig):
    if q.does_file_exist_sync_node(get_save_path(path_s + ".acc.partial")):
        q.qrename_info(get_save_path(path_s + ".acc.partial"), get_save_path(path_s + ".acc"))
    finished_tags = q.properly_truncate_fields(get_save_path(path_s + ".acc"))
    sfw = q.open_fields(get_save_path(path_s + ".acc"), "a", [ 1, 1, 1, 4, ])
    for inv_acc in [ 2, 1 ]:
        for p in wi:
            idx, tslice, inv_type_p, inv_acc_p = p
            if inv_type_p == inv_type and inv_acc_p == inv_acc:
                compute_prop_wsrc(gf, gt, tslice, job_tag, inv_type, inv_acc,
                        idx = idx, sfw = sfw, path_sp = path_sp,
                        psel = psel, fsel = fsel, fselc = fselc, eig = eig,
                        finished_tags = finished_tags)
        q.clean_cache(q.cache_inv)
    sfw.close()
    q.qtouch_info(get_save_path(os.path.join(path_sp, "checkpoint.txt")))
    q.qtouch_info(get_save_path(os.path.join(path_sp, "checkpoint ; wsnk.txt")))
    q.qrename_info(get_save_path(path_s + ".acc"), get_save_path(path_s))
    # q.qar_create_info(get_save_path(path_sp + ".qar"), get_save_path(path_sp), is_remove_folder_after = True)
    # q.qar_create_info(get_save_path(path_s + ".qar"), get_save_path(path_s), is_remove_folder_after = True)

@q.timer
def compute_prop_psrc(gf, gt, xg, job_tag, inv_type, inv_acc, *,
        idx, sfw, sfw_hvp = None, path_sp, psel, fsel, fselc, eig, finished_tags):
    tag = f"xg=({xg[0]},{xg[1]},{xg[2]},{xg[3]}) ; type={inv_type} ; accuracy={inv_acc}"
    if tag in finished_tags:
        return None
    q.check_stop()
    q.check_time_limit()
    q.displayln_info(f"compute_prop_psrc: idx={idx} xg={xg}", job_tag, inv_type, inv_acc)
    inv = ru.get_inv(gf, job_tag, inv_type, inv_acc, eig = eig)
    total_site = ru.get_total_site(job_tag)
    geo = q.Geometry(total_site, 1)
    src = q.mk_point_src(geo, xg)
    fn_sp = os.path.join(path_sp, f"{tag}.lat")
    prop = compute_prop(inv, src, tag = tag, sfw = sfw, fn_sp = fn_sp, psel = psel, fsel = fsel, fselc = fselc)
    if sfw_hvp is not None:
        chvp_16 = q.contract_chvp_16(prop, prop)
        chvp_16.save_float_from_double(sfw_hvp, tag)
        sfw_hvp.flush()
    prop_gt = gt * prop
    fn_spw = os.path.join(path_sp, f"{tag} ; wsnk.lat")
    prop_gt.glb_sum_tslice().save(get_save_path(fn_spw))

def get_common_finished_tags(l1, l2):
    len1 = len(l1)
    len2 = len(l2)
    len_min = min(len1, len2)
    l = l1[:len_min]
    assert l == l2[:len_min]
    return l

@q.timer_verbose
def compute_prop_psrc_all(gf, gt, pi, job_tag, inv_type, *,
        path_s, path_hvp = None, path_sp, psel, fsel, fselc, eig):
    if q.does_file_exist_sync_node(get_save_path(path_s + ".acc.partial")):
        q.qrename_info(get_save_path(path_s + ".acc.partial"), get_save_path(path_s + ".acc"))
    finished_prop_tags = q.properly_truncate_fields(get_save_path(path_s + ".acc"))
    finished_tags = finished_prop_tags
    sfw_hvp = None
    if path_hvp is not None:
        finished_hvp_tags = q.properly_truncate_fields(get_save_path(path_hvp + ".acc"))
        finished_tags = get_common_finished_tags(finished_prop_tags, finished_hvp_tags)
        if finished_tags != finished_prop_tags:
            q.truncate_fields(get_save_path(path_s + ".acc"), finished_tags)
        if finished_tags != finished_hvp_tags:
            q.truncate_fields(get_save_path(path_hvp + ".acc"), finished_tags)
        sfw_hvp = q.open_fields(get_save_path(path_hvp + ".acc"), "a", [ 1, 1, 1, 4, ])
    sfw = q.open_fields(get_save_path(path_s + ".acc"), "a", [ 1, 1, 1, 8, ])
    for inv_acc in [ 2, 1, 0 ]:
        for p in pi:
            idx, xg, inv_type_p, inv_acc_p = p
            if inv_type_p == inv_type and inv_acc_p == inv_acc:
                compute_prop_psrc(gf, gt, xg, job_tag, inv_type, inv_acc,
                        idx = idx, sfw = sfw, sfw_hvp = sfw_hvp, path_sp = path_sp,
                        psel = psel, fsel = fsel, fselc = fselc, eig = eig,
                        finished_tags = finished_tags)
        q.clean_cache(q.cache_inv)
    if sfw_hvp is not None:
        sfw_hvp.close()
    sfw.close()
    q.qtouch_info(get_save_path(os.path.join(path_sp, "checkpoint.txt")))
    q.qrename_info(get_save_path(path_hvp + ".acc"), get_save_path(path_hvp))
    q.qrename_info(get_save_path(path_s + ".acc"), get_save_path(path_s))
    # q.qar_create_info(get_save_path(path_sp + ".qar"), get_save_path(path_sp), is_remove_folder_after = True)
    # q.qar_create_info(get_save_path(path_s + ".qar"), get_save_path(path_s), is_remove_folder_after = True)
    # q.qar_create_info(get_save_path(path_hvp + ".qar"), get_save_path(path_hvp), is_remove_folder_after = True)

@q.timer
def run_prop_psrc_light(job_tag, traj, get_gf, get_eig, get_gt, get_psel, get_fsel, get_pi):
    if None in [ get_gf, get_eig, get_gt, get_psel, get_fsel, get_pi, ]:
        return
    if get_load_path(f"{job_tag}/prop-psrc-light/traj-{traj}/geon-info.txt") is not None:
        return
    if q.obtain_lock(f"locks/{job_tag}-{traj}-psrc-light"):
        gf = get_gf()
        gt = get_gt()
        eig = get_eig()
        fsel, fselc = get_fsel()
        pi = get_pi()
        compute_prop_psrc_all(gf, gt, pi, job_tag, inv_type = 0,
                path_s = f"{job_tag}/prop-psrc-light/traj-{traj}",
                path_hvp = f"{job_tag}/hvp-psrc-light/traj-{traj}",
                path_sp = f"{job_tag}/psel-prop-psrc-light/traj-{traj}",
                psel = get_psel(), fsel = fsel, fselc = fselc, eig = eig)
        q.release_lock()

@q.timer
def run_prop_psrc_strange(job_tag, traj, get_gf, get_gt, get_psel, get_fsel, get_pi):
    if None in [ get_gf, get_gt, get_psel, get_fsel, get_pi, ]:
        return
    if get_load_path(f"{job_tag}/prop-psrc-strange/traj-{traj}/geon-info.txt") is not None:
        return
    if q.obtain_lock(f"locks/{job_tag}-{traj}-psrc-strange"):
        gf = get_gf()
        gt = get_gt()
        fsel, fselc = get_fsel()
        pi = get_pi()
        compute_prop_psrc_all(gf, gt, pi, job_tag, inv_type = 1,
                path_s = f"{job_tag}/prop-psrc-strange/traj-{traj}",
                path_hvp = f"{job_tag}/hvp-psrc-strange/traj-{traj}",
                path_sp = f"{job_tag}/psel-prop-psrc-strange/traj-{traj}",
                psel = get_psel(), fsel = fsel, fselc = fselc, eig = None)
        q.release_lock()

@q.timer
def run_prop_wsrc_light(job_tag, traj, get_gf, get_eig, get_gt, get_psel, get_fsel, get_wi):
    if None in [ get_gf, get_eig, get_gt, get_psel, get_fsel, ]:
        return
    if get_load_path(f"{job_tag}/prop-wsrc-light/traj-{traj}/geon-info.txt") is not None:
        return
    if q.obtain_lock(f"locks/{job_tag}-{traj}-wsrc-light"):
        gf = get_gf()
        gt = get_gt()
        eig = get_eig()
        fsel, fselc = get_fsel()
        wi = get_wi()
        compute_prop_wsrc_all(gf, gt, wi, job_tag, inv_type = 0,
                path_s = f"{job_tag}/prop-wsrc-light/traj-{traj}",
                path_sp = f"{job_tag}/psel-prop-wsrc-light/traj-{traj}",
                psel = get_psel(), fsel = fsel, fselc = fselc, eig = eig)
        q.release_lock()

@q.timer
def run_prop_wsrc_strange(job_tag, traj, get_gf, get_gt, get_psel, get_fsel, get_wi):
    if None in [ get_gf, get_gt, get_psel, get_fsel, ]:
        return
    if get_load_path(f"{job_tag}/prop-wsrc-strange/traj-{traj}/geon-info.txt") is not None:
        return
    if q.obtain_lock(f"locks/{job_tag}-{traj}-wsrc-strange"):
        gf = get_gf()
        gt = get_gt()
        fsel, fselc = get_fsel()
        wi = get_wi()
        compute_prop_wsrc_all(gf, gt, wi, job_tag, inv_type = 1,
                path_s = f"{job_tag}/prop-wsrc-strange/traj-{traj}",
                path_sp = f"{job_tag}/psel-prop-wsrc-strange/traj-{traj}",
                psel = get_psel(), fsel = fsel, fselc = fselc, eig = None)
        q.release_lock()

@q.timer
def run_job(job_tag, traj):
    fns_produce = [
            f"{job_tag}/point-selection/traj-{traj}.txt",
            f"{job_tag}/field-selection/traj-{traj}.field",
            f"{job_tag}/wall-src-info-light/traj-{traj}.txt",
            f"{job_tag}/wall-src-info-strange/traj-{traj}.txt",
            f"{job_tag}/prop-wsrc-strange/traj-{traj}/geon-info.txt",
            f"{job_tag}/prop-wsrc-light/traj-{traj}/geon-info.txt",
            f"{job_tag}/prop-psrc-strange/traj-{traj}/geon-info.txt",
            f"{job_tag}/prop-psrc-light/traj-{traj}/geon-info.txt",
            ]
    fns_need = [
            (f"{job_tag}/configs/ckpoint_lat.{traj}", f"{job_tag}/configs/ckpoint_lat.IEEE64BIG.{traj}",),
            f"{job_tag}/gauge-transform/traj-{traj}.field",
            f"{job_tag}/eig/traj-{traj}",
            f"{job_tag}/eig/traj-{traj}/metadata.txt",
            f"{job_tag}/eig/traj-{traj}/eigen-values.txt",
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
    get_wi = run_wi(job_tag, traj)
    get_pi = run_pi(job_tag, traj, get_psel)
    #
    def run_with_eig():
        get_eig = run_eig(job_tag, traj_gf, get_gf)
        run_get_inverter(job_tag, traj, inv_type = 0, get_gf = get_gf, get_eig = get_eig)
        run_get_inverter(job_tag, traj, inv_type = 0, get_gf = get_gf, get_gt = get_gt, get_eig = get_eig)
        run_prop_wsrc_light(job_tag, traj, get_gf, get_eig, get_gt, get_psel, get_fsel, get_wi)
        run_prop_psrc_light(job_tag, traj, get_gf, get_eig, get_gt, get_psel, get_fsel, get_pi)
    #
    def run_with_eig_strange():
        get_eig = None
        run_get_inverter(job_tag, traj, inv_type = 1, get_gf = get_gf, get_eig = get_eig)
        run_get_inverter(job_tag, traj, inv_type = 1, get_gf = get_gf, get_gt = get_gt, get_eig = get_eig)
        run_prop_wsrc_strange(job_tag, traj, get_gf, get_gt, get_psel, get_fsel, get_wi)
        run_prop_psrc_strange(job_tag, traj, get_gf, get_gt, get_psel, get_fsel, get_pi)
    #
    run_with_eig()
    run_with_eig_strange()
    #
    q.clean_cache()
    q.timer_display()

tag = "trajs"
rup.dict_params["test-4nt8"][tag] = list(range(1000, 1400, 100))
rup.dict_params["test-4nt16"][tag] = list(range(1000, 1400, 100))
rup.dict_params["48I"][tag] = list(range(3000, 500, -5))
rup.dict_params["24D"][tag] = list(range(1000, 10000, 10))
rup.dict_params["16IH2"][tag] = list(range(1000, 10000, 50))
rup.dict_params["32IfineH"][tag] = list(range(1000, 10000, 50))

tag = "n_points_psel"
rup.dict_params["test-4nt8"][tag] = 6
rup.dict_params["test-4nt16"][tag] = 32
rup.dict_params["48I"][tag] = 2048
rup.dict_params["24D"][tag] = 1024
rup.dict_params["32IfineH"][tag] = 512

rup.dict_params["test-4nt8"]["n_points"] = [
        [ 6, 2, 1, ],
        [ 3, 2, 1, ],
        ]

rup.dict_params["test-4nt16"]["n_points"] = [
        [ 32, 4, 2, ],
        [ 16, 4, 2, ],
        ]

rup.dict_params["48I"]["n_points"] = [
        [ 2048, 64, 16, ],
        [ 1024, 64, 16, ],
        ]

# rup.dict_params["test-4nt8"]["fermion_params"][0][2]["Ls"] = 10
# rup.dict_params["test-4nt8"]["fermion_params"][1][2]["Ls"] = 10
# rup.dict_params["test-4nt8"]["fermion_params"][2][2]["Ls"] = 10

# rup.dict_params["test-4nt16"]["fermion_params"][0][2]["Ls"] = 10
# rup.dict_params["test-4nt16"]["fermion_params"][1][2]["Ls"] = 10
# rup.dict_params["test-4nt16"]["fermion_params"][2][2]["Ls"] = 10

tag = "n_exact_wsrc"
rup.dict_params["test-4nt8"][tag] = 2
rup.dict_params["48I"][tag] = 2

tag = "prob_exact_wsrc"
rup.dict_params["test-4nt16"][tag] = 1/8
rup.dict_params["16IH2"][tag] = 1/16
rup.dict_params["32IfineH"][tag] = 1/32

qg.begin_with_gpt()

# ADJUST ME
job_tags = [
        "test-4nt8", "test-4nt16",
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
