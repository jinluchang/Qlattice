#!/usr/bin/env python3

import qlat as q
import rbc_ukqcd as ru
import rbc_ukqcd_params as rup
import pprint

import os

def get_save_path(fn):
    return os.path.join("results", fn)

def get_load_path(fn):
    if fn is None:
        return None
    path_list = [
            "results",
            "../mk-gf-gt/results",
            "../mk-lanc/results",
            "../mk-selected-data/results",
            "/sdcc/u/jluchang/qcdqedta/hlbl-data",
            ]
    for path in path_list:
        p = os.path.join(path, fn)
        if q.does_file_exist_sync_node(p):
            return p
    return None

@q.timer_verbose
def check_job(job_tag, traj):
    # return True if config is finished or unavailable
    fns_produce = [
            get_load_path(f"psel-prop-wsrc-light/{job_tag}/traj={traj}/checkpoint ; wsnk.txt"),
            get_load_path(f"psel-prop-wsrc-strange/{job_tag}/traj={traj}/checkpoint ; wsnk.txt"),
            ]
    is_job_done = True
    for fn in fns_produce:
        if fn is None:
            q.displayln_info(f"check_job: {job_tag} {traj} to do as some file does not exist.")
            is_job_done = False
            break
    if is_job_done:
        return True
    #
    fns_need = [
            get_load_path(f"point-selection/{job_tag}/traj={traj}.txt"),
            get_load_path(f"field-selection/{job_tag}/traj={traj}.field"),
            get_load_path(f"prop-wsrc-strange/{job_tag}/traj={traj}"),
            get_load_path(f"prop-wsrc-light/{job_tag}/traj={traj}"),
            get_load_path(f"gauge-transform/{job_tag}/traj={traj}.field"),
            ]
    for fn in fns_need:
        if fn is None:
            q.displayln_info(f"check_job: {job_tag} {traj} unavailable as {fn} does not exist.")
            return True
    #
    q.check_stop()
    q.check_time_limit()
    #
    q.qmkdir_info(f"locks")
    q.qmkdir_info(get_save_path(f""))
    #
    return False

@q.timer_verbose
def compute_prop_wsrc_all(path_sp, path_s, fsel, fselc):
    q.check_stop()
    q.check_time_limit()
    total_site = fsel.total_site()
    psel = q.get_psel_tslice(total_site)
    prob = fsel.prob()
    sprop_c = q.SelProp(fselc)
    sprop = q.SelProp(fsel)
    sp_prop = q.PselProp(psel)
    sfr = q.open_fields(get_load_path(path_s), "r")
    tags = sfr.list()
    for tag in tags:
        fn_spw = os.path.join(path_sp, f"{tag} ; wsnk.lat")
        if q.does_file_exist_sync_node(fn_spw):
            continue
        sprop_c.load_double_from_float(sfr, tag)
        sprop @= sprop_c
        sp_prop.swap(sprop.glb_sum_tslice())
        sp_prop *= 1 / prob
        sp_prop.save(get_save_path(fn_spw))
    sfr.close()
    q.qtouch(get_save_path(f"psel-prop-wsrc-light/{job_tag}/traj={traj}/checkpoint ; wsnk.txt"))

rup.dict_params["test-4nt8"]["trajs"] = list(range(500, 3000, 5))
rup.dict_params["test-4nt16"]["trajs"] = list(range(500, 3000, 5))
rup.dict_params["24D"]["trajs"] = list(range(500, 3000, 5))
rup.dict_params["48I"]["trajs"] = list(range(500, 3000, 5))
rup.dict_params["64I"]["trajs"] = list(range(500, 3000, 5))

# ADJUST ME
job_tags = [
        "test-4nt8",
        "test-4nt16",
        # "24D",
        # "48I",
        # "64I",
        ]

q.check_time_limit()

for job_tag in job_tags:
    q.displayln_info(pprint.pformat(rup.dict_params[job_tag]))
    for traj in rup.dict_params[job_tag]["trajs"]:
        run_job(job_tag, traj)



