#!/usr/bin/env python3

# Need --mpi X.X.X.X --mpi X.X.X runtime option

import qlat_gpt as qg
import rbc_ukqcd as ru

from jobs import *

load_path_list[:] = [
        "results",
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
def run_prop_wsrc_strange(job_tag, traj, get_gf, get_eig, get_gt, get_psel, get_fsel, get_wi):
    if None in [ get_gf, get_eig, get_gt, get_psel, get_fsel, ]:
        return
    if get_load_path(f"prop-wsrc-strange/{job_tag}/traj={traj}") is not None:
        return
    if q.obtain_lock(f"locks/{job_tag}-{traj}-wsrc-strange"):
        gf = get_gf()
        gt = get_gt()
        eig = get_eig()
        fsel, fselc = get_fsel()
        wi = get_wi()
        compute_prop_wsrc_all(gf, gt, wi, job_tag, inv_type = 1,
                path_s = f"prop-wsrc-strange/{job_tag}/traj={traj}",
                path_sp = f"psel-prop-wsrc-strange/{job_tag}/traj={traj}",
                psel = get_psel(), fsel = fsel, fselc = fselc, eig = eig)
        q.release_lock()

@q.timer_verbose
def run_job(job_tag, traj):
    fns_produce = [
            f"wall-src-info-light/{job_tag}/traj={traj}.txt",
            f"wall-src-info-strange/{job_tag}/traj={traj}.txt",
            f"prop-wsrc-strange/{job_tag}/traj={traj}",
            f"prop-wsrc-light/{job_tag}/traj={traj}",
            ]
    fns_need = [
            (f"configs/{job_tag}/ckpoint_lat.{traj}", f"configs/{job_tag}/ckpoint_lat.IEEE64BIG.{traj}",),
            f"gauge-transform/{job_tag}/traj={traj}.field",
            f"point-selection/{job_tag}/traj={traj}.txt",
            f"field-selection/{job_tag}/traj={traj}.field",
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
    assert get_wi is not None
    #
    def run_with_eig():
        get_eig = run_eig(job_tag, traj_gf, get_gf)
        run_prop_wsrc_light(job_tag, traj, get_gf, get_eig, get_gt, get_psel, get_fsel, get_wi)
    #
    run_with_eig()
    #
    def run_with_eig_strange():
        get_eig_strange = run_eig_strange(job_tag, traj_gf, get_gf)
        run_prop_wsrc_strange(job_tag, traj, get_gf, get_eig_strange, get_gt, get_psel, get_fsel, get_wi)
    #
    run_with_eig_strange()
    #
    q.clean_cache()
    q.timer_display()

rup.dict_params["test-4nt8"]["trajs"] = list(range(1000, 1400, 100))
rup.dict_params["test-4nt16"]["trajs"] = list(range(1000, 1400, 100))
rup.dict_params["48I"]["trajs"] = list(range(3000, 500, -5))
rup.dict_params["16IH2"]["trajs"] = list(range(1000, 10000, 50))
rup.dict_params["32IfineH"]["trajs"] = list(range(1000, 10000, 50))

rup.dict_params["test-4nt8"]["fermion_params"][0][2]["Ls"] = 10
rup.dict_params["test-4nt8"]["fermion_params"][1][2]["Ls"] = 10
rup.dict_params["test-4nt8"]["fermion_params"][2][2]["Ls"] = 10

# rup.dict_params["test-4nt16"]["fermion_params"][0][2]["Ls"] = 10
# rup.dict_params["test-4nt16"]["fermion_params"][1][2]["Ls"] = 10
# rup.dict_params["test-4nt16"]["fermion_params"][2][2]["Ls"] = 10

rup.dict_params["test-4nt8"]["n_exact_wsrc"] = 2
rup.dict_params["48I"]["n_exact_wsrc"] = 2

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
        # "16IH2",
        # "32IfineH",
        # "48I",
        ]

q.check_time_limit()

for job_tag in job_tags:
    q.displayln_info(pprint.pformat(rup.dict_params[job_tag]))
    for traj in rup.dict_params[job_tag]["trajs"]:
        run_job(job_tag, traj)

qg.end_with_gpt()
