#!/usr/bin/env python3

# Need --mpi X.X.X.X --mpi X.X.X runtime option

import qlat_gpt as qg
import rbc_ukqcd as ru

from jobs import *

load_path_list[:] = [
        "results",
        "/gpfs/alpine/lgt116/proj-shared/ljin",
        "../mk-gf-gt/results",
        "../mk-sel/results",
        "../mk-lanc/results",
        "../qcddata",
        os.path.join(os.getenv("HOME"), "Qlat-sample-data/mk-gf-gt/results"),
        os.path.join(os.getenv("HOME"), "Qlat-sample-data/mk-sel/results"),
        os.path.join(os.getenv("HOME"), "Qlat-sample-data/mk-lanc/results"),
        os.path.join(os.getenv("HOME"), "qcddata"),
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

@q.timer
def run_prop_psrc_light(job_tag, traj, get_gf, get_eig, get_gt, get_psel, get_fsel):
    if None in [ get_gf, get_eig, get_gt, get_psel, get_fsel, ]:
        return
    if get_load_path(f"prop-psrc-light/{job_tag}/traj={traj}") is not None:
        return
    if q.obtain_lock(f"locks/{job_tag}-{traj}-psrc-light"):
        gf = get_gf()
        gt = get_gt()
        eig = get_eig()
        fsel, fselc = get_fsel()
        pi = mk_rand_point_src_info(job_tag, traj, get_psel())
        save_point_src_info(pi, get_save_path(f"point-src-info/{job_tag}/traj={traj}.txt"));
        compute_prop_psrc_all(gf, gt, pi, job_tag, inv_type = 0,
                path_s = f"prop-psrc-light/{job_tag}/traj={traj}",
                path_hvp = f"hvp-psrc-light/{job_tag}/traj={traj}",
                path_sp = f"psel-prop-psrc-light/{job_tag}/traj={traj}",
                psel = get_psel(), fsel = fsel, fselc = fselc, eig = eig)
        q.release_lock()

@q.timer
def run_prop_psrc_strange(job_tag, traj, get_gf, get_gt, get_psel, get_fsel):
    if None in [ get_gf, get_gt, get_psel, get_fsel, ]:
        return
    if get_load_path(f"prop-psrc-strange/{job_tag}/traj={traj}") is not None:
        return
    if q.obtain_lock(f"locks/{job_tag}-{traj}-psrc-strange"):
        gf = get_gf()
        gt = get_gt()
        fsel, fselc = get_fsel()
        pi = mk_rand_point_src_info(job_tag, traj, get_psel())
        save_point_src_info(pi, get_save_path(f"point-src-info/{job_tag}/traj={traj}.txt"));
        compute_prop_psrc_all(gf, gt, pi, job_tag, inv_type = 1,
                path_s = f"prop-psrc-strange/{job_tag}/traj={traj}",
                path_hvp = f"hvp-psrc-strange/{job_tag}/traj={traj}",
                path_sp = f"psel-prop-psrc-strange/{job_tag}/traj={traj}",
                psel = get_psel(), fsel = fsel, fselc = fselc, eig = None)
        q.release_lock()

@q.timer
def run_prop_wsrc_light(job_tag, traj, get_gf, get_eig, get_gt, get_psel, get_fsel, get_wi):
    if None in [ get_gf, get_eig, get_gt, get_psel, get_fsel, ]:
        return
    if get_load_path(f"prop-wsrc-light/{job_tag}/traj={traj}") is not None:
        return
    if q.obtain_lock(f"locks/{job_tag}-{traj}-wsrc-light"):
        gf = get_gf()
        gt = get_gt()
        eig = get_eig()
        fsel, fselc = get_fsel()
        wi = get_wi()
        compute_prop_wsrc_all(gf, gt, wi, job_tag, inv_type = 0,
                path_s = f"prop-wsrc-light/{job_tag}/traj={traj}",
                path_sp = f"psel-prop-wsrc-light/{job_tag}/traj={traj}",
                psel = get_psel(), fsel = fsel, fselc = fselc, eig = eig)
        q.release_lock()

@q.timer
def run_prop_wsrc_strange(job_tag, traj, get_gf, get_gt, get_psel, get_fsel, get_wi):
    if None in [ get_gf, get_gt, get_psel, get_fsel, ]:
        return
    if get_load_path(f"prop-wsrc-strange/{job_tag}/traj={traj}") is not None:
        return
    if q.obtain_lock(f"locks/{job_tag}-{traj}-wsrc-strange"):
        gf = get_gf()
        gt = get_gt()
        fsel, fselc = get_fsel()
        wi = get_wi()
        compute_prop_wsrc_all(gf, gt, wi, job_tag, inv_type = 1,
                path_s = f"prop-wsrc-strange/{job_tag}/traj={traj}",
                path_sp = f"psel-prop-wsrc-strange/{job_tag}/traj={traj}",
                psel = get_psel(), fsel = fsel, fselc = fselc, eig = None)
        q.release_lock()

@q.timer
def run_job(job_tag, traj):
    fns_produce = [
            f"point-selection/{job_tag}/traj={traj}.txt",
            f"field-selection/{job_tag}/traj={traj}.field",
            f"wall-src-info-light/{job_tag}/traj={traj}.txt",
            f"wall-src-info-strange/{job_tag}/traj={traj}.txt",
            f"prop-wsrc-strange/{job_tag}/traj={traj}",
            f"prop-wsrc-light/{job_tag}/traj={traj}",
            f"prop-psrc-strange/{job_tag}/traj={traj}",
            f"prop-psrc-light/{job_tag}/traj={traj}",
            ]
    fns_need = [
            (f"configs/{job_tag}/ckpoint_lat.{traj}", f"configs/{job_tag}/ckpoint_lat.IEEE64BIG.{traj}",),
            f"gauge-transform/{job_tag}/traj={traj}.field",
            f"eig/{job_tag}/traj={traj}",
            f"eig/{job_tag}/traj={traj}/metadata.txt",
            f"eig/{job_tag}/traj={traj}/eigen-values.txt",
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
    #
    def run_with_eig():
        get_eig = run_eig(job_tag, traj_gf, get_gf)
        run_prop_wsrc_light(job_tag, traj, get_gf, get_eig, get_gt, get_psel, get_fsel, get_wi)
        run_prop_psrc_light(job_tag, traj, get_gf, get_eig, get_gt, get_psel, get_fsel)
    #
    run_with_eig()
    #
    run_prop_wsrc_strange(job_tag, traj, get_gf, get_gt, get_psel, get_fsel, get_wi)
    run_prop_psrc_strange(job_tag, traj, get_gf, get_gt, get_psel, get_fsel)
    #
    q.clean_cache()
    q.timer_display()

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

rup.dict_params["test-4nt8"]["trajs"] = list(range(1000, 1400, 100))
rup.dict_params["test-4nt16"]["trajs"] = list(range(1000, 1400, 100))
rup.dict_params["48I"]["trajs"] = list(range(3000, 500, -5))

rup.dict_params["test-4nt8"]["fermion_params"][0][2]["Ls"] = 10
rup.dict_params["test-4nt8"]["fermion_params"][1][2]["Ls"] = 10
rup.dict_params["test-4nt8"]["fermion_params"][2][2]["Ls"] = 10

# rup.dict_params["test-4nt16"]["fermion_params"][0][2]["Ls"] = 10
# rup.dict_params["test-4nt16"]["fermion_params"][1][2]["Ls"] = 10
# rup.dict_params["test-4nt16"]["fermion_params"][2][2]["Ls"] = 10

rup.dict_params["test-4nt8"]["n_exact_wsrc"] = 2
rup.dict_params["48I"]["n_exact_wsrc"] = 2

rup.dict_params["test-4nt16"]["prob_exact_wsrc"] = 1/8
rup.dict_params["16IH2"]["prob_exact_wsrc"] = 1/16
rup.dict_params["32IfineH"]["prob_exact_wsrc"] = 1/32

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
