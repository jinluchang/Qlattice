#!/usr/bin/env python3

# Need --mpi X.X.X.X --mpi X.X.X runtime option

import qlat_grid_io as qgi

from jobs import *

load_path_list[:] = [
        "results",
        "../qcddata",
        os.path.join(os.getenv("HOME"), "qcddata"),
        "../mk-gf-gt/results",
        "../mk-sel/results",
        "../mk-lanc/results",
        os.path.join(os.getenv("HOME"), "Qlat-sample-data/default/mk-gf-gt/results"),
        os.path.join(os.getenv("HOME"), "Qlat-sample-data/default/mk-sel/results"),
        os.path.join(os.getenv("HOME"), "Qlat-sample-data/default/mk-lanc/results"),
        os.path.join(os.getenv("HOME"), "Qlat-sample-data/default/mk-wsrc-prop/results"),
        "/sdcc/u/jluchang/qcdqedta/summit-oakforest-data-cache",
        ]

@q.timer_verbose
def run_prop_wsrc_light(job_tag, traj, get_gf, get_gt, get_psel, get_fsel, get_wi):
    if None in [ get_gf, get_gt, get_psel, get_fsel, ]:
        return
    path = get_load_path(f"{job_tag}/prop-wsrc-light/traj-{traj}/geon-info.txt")
    if path is None:
        q.displayln_info(f"run_prop_wsrc_light: {job_tag} {traj} props do not exist")
        return
    if q.obtain_lock(f"locks/{job_tag}-{traj}-wsrc-light"):
        gt = get_gt()
        geo = gt.geo()
        gt_inv = gt.inv()
        fsel, fselc = get_fsel()
        wi = get_wi()
        sfr = q.open_fields(path, "r")
        tag_list = sfr.list()
        for p in wi:
            idx, tslice, inv_type, inv_acc = p
            fn = os.path.join(f"{job_tag}/prop-wsrc-light-scidac-full/traj-{traj}",
                    f"tslice-{tslice}_type-{inv_type}_accuracy-{inv_acc}.prop.scidac")
            if get_load_path(fn) is not None:
                q.displayln_info(f"'{get_load_path(fn)}' exist")
                continue
            tag = f"tslice={tslice} ; type={inv_type} ; accuracy={inv_acc}"
            if tag not in tag_list:
                q.displayln_info(f"run_prop_wsrc_light: {job_tag} {traj} tag={tag} not in data")
            sc_prop = q.SelProp(fselc)
            sc_prop.load_double_from_float(sfr, tag)
            sc_prop = gt_inv * sc_prop
            s_prop = q.SelProp(fsel)
            s_prop @= sc_prop
            prop = q.Prop(geo)
            q.set_zero(prop)
            prop @= s_prop
            qgi.save_prop_float(prop, get_save_path(fn))
        sfr.close()
        q.qtouch_info(get_save_path(f"{job_tag}/prop-wsrc-light-scidac-full/traj-{traj}/checkpoint.txt"))
        q.release_lock()

@q.timer_verbose
def run_prop_wsrc_strange(job_tag, traj, get_gf, get_gt, get_psel, get_fsel, get_wi):
    if None in [ get_gf, get_gt, get_psel, get_fsel, ]:
        return
    path = get_load_path(f"{job_tag}/prop-wsrc-strange/traj-{traj}/geon-info.txt")
    if path is None:
        q.displayln_info(f"run_prop_wsrc_strange: {job_tag} {traj} props do not exist")
        return
    if q.obtain_lock(f"locks/{job_tag}-{traj}-wsrc-strange"):
        gt = get_gt()
        geo = gt.geo()
        gt_inv = gt.inv()
        fsel, fselc = get_fsel()
        wi = get_wi()
        sfr = q.open_fields(path, "r")
        tag_list = sfr.list()
        for p in wi:
            idx, tslice, inv_type, inv_acc = p
            fn = os.path.join(f"{job_tag}/prop-wsrc-strange-scidac-full/traj-{traj}",
                    f"tslice-{tslice}_type-{inv_type}_accuracy-{inv_acc}.prop.scidac")
            if get_load_path(fn) is not None:
                q.displayln_info(f"'{get_load_path(fn)}' exist")
                continue
            tag = f"tslice={tslice} ; type={inv_type} ; accuracy={inv_acc}"
            if tag not in tag_list:
                q.displayln_info(f"run_prop_wsrc_strange: {job_tag} {traj} tag={tag} not in data")
            sc_prop = q.SelProp(fselc)
            sc_prop.load_double_from_float(sfr, tag)
            sc_prop = gt_inv * sc_prop
            s_prop = q.SelProp(fsel)
            s_prop @= sc_prop
            prop = q.Prop(geo)
            q.set_zero(prop)
            prop @= s_prop
            qgi.save_prop_float(prop, get_save_path(fn))
        sfr.close()
        q.qtouch_info(get_save_path(f"{job_tag}/prop-wsrc-strange-scidac-full/traj-{traj}/checkpoint.txt"))
        q.release_lock()

@q.timer_verbose
def run_job(job_tag, traj):
    fns_produce = [
            f"{job_tag}/prop-wsrc-light-scidac-full/traj-{traj}/checkpoint.txt",
            f"{job_tag}/prop-wsrc-strange-scidac-full/traj-{traj}/checkpoint.txt",
            ]
    fns_need = [
            # (f"{job_tag}/configs/ckpoint_lat.{traj}", f"{job_tag}/configs/ckpoint_lat.IEEE64BIG.{traj}",),
            f"{job_tag}/gauge-transform/traj-{traj}.field",
            f"{job_tag}/point-selection/traj-{traj}.txt",
            f"{job_tag}/field-selection/traj-{traj}.field",
            f"{job_tag}/wall-src-info-light/traj-{traj}.txt",
            f"{job_tag}/wall-src-info-strange/traj-{traj}.txt",
            f"{job_tag}/prop-wsrc-strange/traj-{traj}/geon-info.txt",
            f"{job_tag}/psel-prop-wsrc-strange/traj-{traj}/checkpoint.txt",
            f"{job_tag}/prop-wsrc-light/traj-{traj}/geon-info.txt",
            f"{job_tag}/psel-prop-wsrc-light/traj-{traj}/checkpoint.txt",
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
    run_prop_wsrc_light(job_tag, traj, get_gf, get_gt, get_psel, get_fsel, get_wi)
    run_prop_wsrc_strange(job_tag, traj, get_gf, get_gt, get_psel, get_fsel, get_wi)
    #
    q.clean_cache()
    q.timer_display()

tag = "trajs"
rup.dict_params["test-4nt8"][tag] = list(range(1000, 1400, 100))
rup.dict_params["test-4nt16"][tag] = list(range(1000, 1400, 100))
rup.dict_params["48I"][tag] = list(range(1000, 3000, 5))
rup.dict_params["24D"][tag] = list(range(1000, 10000, 10))
rup.dict_params["24DH"][tag] = list(range(200, 1000, 10))
rup.dict_params["32Dfine"][tag] = list(range(1000, 10000, 10))
rup.dict_params["16IH2"][tag] = list(range(1000, 10000, 10))
rup.dict_params["32IfineH"][tag] = list(range(1000, 10000, 10))
rup.dict_params["32IcoarseH1"][tag] = list(range(300, 2000, 50))
rup.dict_params["24IH1"][tag] = list(range(1000, 10000, 100))
rup.dict_params["24IH2"][tag] = list(range(1000, 10000, 100))
rup.dict_params["24IH3"][tag] = list(range(1000, 10000, 100))
rup.dict_params["24IH4"][tag] = list(range(1000, 10000, 100))
rup.dict_params["32IH1"][tag] = list(range(1000, 10000, 50))
rup.dict_params["32IH2"][tag] = list(range(1000, 10000, 100)) + list(range(1040, 10000, 100))
rup.dict_params["32IH3"][tag] = list(range(1000, 10000, 50))

rup.dict_params["test-4nt8"]["fermion_params"][0][2]["Ls"] = 10
rup.dict_params["test-4nt8"]["fermion_params"][1][2]["Ls"] = 10
rup.dict_params["test-4nt8"]["fermion_params"][2][2]["Ls"] = 10

# rup.dict_params["test-4nt16"]["fermion_params"][0][2]["Ls"] = 10
# rup.dict_params["test-4nt16"]["fermion_params"][1][2]["Ls"] = 10
# rup.dict_params["test-4nt16"]["fermion_params"][2][2]["Ls"] = 10

qgi.begin_with_grid()

# ADJUST ME
job_tags = [
        "test-4nt8", "test-4nt16",
        # "32IH1",
        # "32IH2",
        # "24IH1",
        # "24IH2",
        # "24IH3",
        # "24D",
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

qgi.end_with_grid()
