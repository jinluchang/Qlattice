#!/usr/bin/env python3

# Need --mpi X.X.X.X runtime option

import qlat as q
import gpt as g
import qlat_gpt as qg
import rbc_ukqcd as ru

import os

def get_save_path(fn):
    return os.path.join("results", fn)

def get_load_path(fn):
    if fn is None:
        return None
    path_list = [ "results",
            "../mk-sample-selected-data/results",
            "/home/luchang/application/Public/Qlat-Grid-cc/sample-data/mk-sample-selected-data/results" ]
    for path in path_list:
        p = os.path.join(path, fn)
        if q.does_file_exist_sync_node(p):
            return p
    return None

@q.timer
def load_point_src_info(path):
    # pi is a list of [ idx xg inv_type inv_acc ]
    dt = q.qload_datatable_sync_node(path, True)
    t = [ list(map(int, l)) for l in dt ][1:]
    pi = [ [ l[0], l[1:5], l[5], l[6] ] for l in t ]
    return pi

@q.timer
def compute_light_hvp(job_tag, traj, *, pi):
    pass

@q.timer
def check_job(job_tag, traj):
    # return True if config is finished
    fns_require = []
    fns_require.append(get_load_path(f"configs/{job_tag}/ckpoint_lat.{traj}"))
    fns_require.append(get_load_path(f"point-src-info/{job_tag}/traj={traj}.txt"))
    for fn in fns_require:
        if fn is None:
            return True
    fns = []
    fns.append(get_load_path(f"hvp/{job_tag}/traj={traj}.lat"))
    for fn in fns:
        if fn is None:
            return False
    return True

@q.timer
def run_job(job_tag, traj):
    if check_job(job_tag, traj):
        return
    q.qmkdir_info("locks")
    q.qmkdir_info(get_save_path(f""))
    q.qmkdir_info(get_save_path(f"hvp"))
    q.qmkdir_info(get_save_path(f"hvp/{job_tag}"))
    #
    total_site = ru.get_total_site(job_tag)
    geo = q.Geometry(total_site, 1)
    q.displayln_info("geo.show() =", geo.show())
    #
    path_gf = get_load_path(f"configs/{job_tag}/ckpoint_lat.{traj}")
    gf = q.GaugeField()
    gf.load(path_gf)
    #
    pi = load_point_src_info(get_load_path(f"point-src-info/{job_tag}/traj={traj}.txt"))
    for p in pi:
        q.displayln_info(p)

qg.begin_with_gpt()

for job_tag in [ "test-4nt16" ]:
    for traj in range(1000, 1100, 100):
        run_job(job_tag, traj)
        q.timer_display()

qg.end_with_gpt()
