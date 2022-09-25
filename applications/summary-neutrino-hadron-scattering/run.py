#!/usr/bin/env python3

import os
import sys
import qlat as q
import numpy as np
import h5py

@q.timer
def get_load_path(fn):
    if fn is None:
        return None
    path_list = [
            "data",
            "../contract-wall-point/analysis",
            "/home/luchang/application/Public/Qlat-Grid-cc/sample-data/contract-wall-point/analysis",
            "/home/frank/application/Public/Qlat-Grid-cc/sample-data/contract-wall-point/analysis",
            "/hpcgpfs01/work/lqcd/qcdqedta/luchang/data-analysis/analysis",
            ]
    for path in path_list:
        if fn == "":
            p = path
        else:
            p = os.path.join(path, fn)
        if q.does_file_exist_sync_node(p):
            return p
    return None

def get_save_path(fn):
    return os.path.join("results", fn)

@q.timer
def check_job(job_tag, traj):
    if q.does_file_exist_sync_node(get_save_path(f"{job_tag}/traj-{traj}.hdf5")):
        q.displayln_info(f"check_job: Job finished {job_tag} {traj}")
        return False
    if get_load_path(f"{job_tag}/field-meson-vv-meson/results-{traj}/checkpoint.txt") is None:
        return False
    if get_load_path(f"{job_tag}/lat-three-point/results-{traj}/checkpoint.txt") is None:
        return False
    if get_load_path(f"{job_tag}/lat-two-point/results-{traj}/checkpoint.txt") is None:
        return False
    return True

def get_t_size(job_tag):
    if job_tag in [ "24D", "32D", "24DH", "32Dfine", ]:
        return 64
    elif job_tag == "48I":
        return 96
    elif job_tag == "64I":
        return 128
    elif job_tag == "test-4nt16":
        return 16
    else:
        assert False

def get_tsep(job_tag):
    if job_tag in [ "24D", "32D", "24DH", ]:
        return 8
    elif job_tag == "32Dfine":
        return 10
    elif job_tag == "48I":
        return 12
    elif job_tag == "64I":
        return 18
    elif job_tag == "test-4nt16":
        return 2
    else:
        assert False

def get_analysis_lmom_list():
    return [
            [ 0, 0, 0, 0, ],
            [ 1, 0, 0, 0, ],
            [-1, 0, 0, 0, ],
            [ 0, 1, 0, 0, ],
            [ 0,-1, 0, 0, ],
            [ 0, 0, 1, 0, ],
            [ 0, 0,-1, 0, ],
            [ 2, 0, 0, 0, ],
            [-2, 0, 0, 0, ],
            [ 0, 2, 0, 0, ],
            [ 0,-2, 0, 0, ],
            [ 0, 0, 2, 0, ],
            [ 0, 0,-2, 0, ],
            ]

@q.timer_verbose
def analysis_two_point(job_tag, traj, type1, type2):
    ld = q.LatData()
    ld.load(get_load_path(f"{job_tag}/lat-two-point/results-{traj}/two-point-wall-snk-sparse-corrected-{type1}-{type2}.lat"))
    arr = np.array([ ld[[t, 15, 15,]][0] for t in range(get_t_size(job_tag)) ])
    # q.displayln_info(arr)
    return arr

@q.timer_verbose
def analysis_three_point(job_tag, traj, type1, type2, type3):
    ld = q.LatData()
    ld.load(get_load_path(f"{job_tag}/lat-three-point/results-{traj}/three-point-{type1}-{type2}-{type3}.lat"))
    arr = np.array([
        [ ld[[t, top, 8,]][0] for top in range(get_t_size(job_tag)) ]
        for t in range(get_t_size(job_tag)) ])
    # q.displayln_info(arr)
    return arr

@q.timer_verbose
def analysis_meson_vv_meson(job_tag, traj, type1, type2, type3, type4):
    f = q.Field("Complex")
    f.load_double_from_float(
            get_load_path(
                f"{job_tag}/field-meson-vv-meson/results-{traj}/forward-{type1}-{type2}-{type3}-{type4}.field"))
    geo = f.geo()
    data = []
    for lmom in get_analysis_lmom_list():
        lmom_factor = list(map(lambda x: -x, lmom))
        f_factor = q.mk_phase_field(geo, lmom_factor)
        f1 = f.copy()
        f1 *= f_factor
        corr = f1.glb_sum_tslice()
        data.append(corr)
    data = np.array(data)
    n_lmom, t_size, n_elem = data.shape
    assert n_lmom == len(get_analysis_lmom_list())
    assert t_size == get_t_size(job_tag)
    assert n_elem == 64
    arr = data.reshape((n_lmom, t_size, 8, 8,))
    # q.displayln_info(arr[0][2])
    return arr

@q.timer_verbose
def analysis(job_tag, traj):
    q.check_stop()
    q.check_time_limit()
    if q.obtain_lock(f"locks/{job_tag}-{traj}"):
        q.timer_reset()
        q.displayln_info(f"analysis {job_tag} {traj}")
        q.qmkdir_info(get_save_path(f"{job_tag}"))
        data = {}
        data["lmom-list"] = np.array(get_analysis_lmom_list())
        data["two-point-0-0"] = analysis_two_point(job_tag, traj, 0, 0)
        data["three-point-0-0-0"] = analysis_three_point(job_tag, traj, 0, 0, 0)
        data["field-meson-vv-meson-0-0-0-0"] = analysis_meson_vv_meson(job_tag, traj, 0, 0, 0, 0)
        if q.get_id_node() == 0:
            with h5py.File(get_save_path(f"{job_tag}/traj-{traj}.hdf5"), "w") as f:
                for k, v in data.items():
                    f.create_dataset(k, data = v)
        q.sync_node()
        q.timer_display()
        q.release_lock()

job_tags = [
        "24D",
        "48I",
        "64I",
        "32D",
        "32Dfine",
        "24DH",
        "test-4nt16",
        ]

trajs = list(range(200, 3000, 10))

size_node_list = [
    [1, 1, 1, 1],
    [1, 1, 1, 2],
    [1, 1, 1, 4],
    [1, 1, 1, 8],
    [1, 1, 1, 16],
    ]

q.begin(sys.argv, size_node_list)

q.qmkdir_info("locks")
q.qmkdir_info("results")

for job_tag in job_tags:
    for traj in trajs:
        if check_job(job_tag, traj):
            analysis(job_tag, traj)

q.timer_display()

q.end()
