#!/usr/bin/env python3

import os
import sys
import qlat_ext as q

@q.timer
def get_load_path(fn):
    if fn is None:
        return None
    path_list = [
            "data",
            "../mk-sample-selected-data/results",
            "/home/luchang/application/Public/Qlat-Grid-cc/sample-data/mk-sample-selected-data/results",
            "/hpcgpfs01/work/lqcd/qcdqedta/hlbl-data-with-cache",
            ]
    for path in path_list:
        if fn == "":
            p = path
        else:
            p = os.path.join(path, fn)
        if q.does_file_exist_sync_node(p):
            return p
    return None

@q.timer
def check_job(job_tag, traj):
    if get_load_path(f"gauge-transform/{job_tag}/traj={traj}.field") is None:
        return False
    if get_load_path(f"prop-psrc-light/{job_tag}/traj={traj}") is None:
        return False
    if get_load_path(f"prop-psrc-strange/{job_tag}/traj={traj}") is None:
        return False
    if get_load_path(f"prop-wsrc-light/{job_tag}/traj={traj}") is None:
        return False
    if get_load_path(f"prop-wsrc-strange/{job_tag}/traj={traj}") is None:
        return False
    return True

@q.timer_verbose
def compute(job_tag, traj):
    q.setup(job_tag)
    q.compute_meson_vv(job_tag, traj)
    q.compute_meson_vv_meson(job_tag, traj)
    q.compute_meson_snk_src(job_tag, traj)
    q.compute_meson_chvp(job_tag, traj)
    q.compute_chvp(job_tag, traj)
    q.compute_two_point_func(job_tag, traj)
    q.compute_three_point_func(job_tag, traj)
    q.compute_psel_fsel_distribution(job_tag, traj)
    q.clear_all_data_cache()
    q.timer_display()

job_tags = [
        "test-4nt16",
        "24D",
        "24DH",
        "32D",
        "32Dfine",
        ]
trajs = list(range(200, 3000, 10))

size_node_list = [
    [1, 1, 1, 1],
    [1, 1, 1, 2],
    [1, 1, 1, 4],
    [1, 1, 1, 8],
    [2, 2, 2, 2],
    [2, 2, 2, 4]]

q.begin(sys.argv, size_node_list)

q.setup_log_idx()
q.setup()

q.qmkdir_info("locks")
q.qmkdir_info("analysis")

q.set_data_path(get_load_path(""))

for job_tag in job_tags:
    for traj in trajs:
        if check_job(job_tag, traj):
            compute(job_tag, traj)

q.timer_display()

q.end()
