#!/usr/bin/env python3

import sys, os
import numpy as np
import glob
import pickle
import qlat as q

import qlat_scripts.v1 as qs

from qlat_scripts.v1 import (
    set_param,
    get_param,
)

usage = f"""
{__file__} --grid XX.XX.XX.XX --src PATH_SRC --dst PATH_DST
# Convert all "*.txt" point-selection files in ``PATH_SRC'' to "*.lati" files in ``PATH_DST'' with ``total_site'' given by ``--grid''.
# E.g.: {__file__} --grid 4.4.4.8 --src results/test-4nt8/point-selection --dst results/test-4nt8/points-selection
#
{__file__}
# Generate some test data and then perform the conversion.
#
"""

def parse_coordinate(x_str):
    x_str_list = x_str.split(".")
    x_list = [ int(s) for s in x_str_list ]
    x = q.Coordinate(x_list)
    return x

@q.timer
def mk_psel(total_site, rate, rs):
    geo = q.Geometry(total_site)
    f_weight = q.FieldRealD(geo, 1)
    f_weight.set_unit()
    f_rand_01 = q.FieldRealD(geo, 1)
    f_rand_01.set_rand(rs, 1.0, 0.0)
    psel = qs.make_psel_from_weight(f_weight, f_rand_01, rate)
    return psel

@q.timer(is_timer_fork=True)
def gen_test_data():
    total_site = q.Coordinate([ 4, 4, 4, 8, ])
    path_src = "results/test-4nt8/point-selection"
    path_dst = "results/test-4nt8/points-selection"
    rate = 0.2
    traj_list = list(range(1000, 2000, 100))
    for traj in traj_list:
        rs = q.RngState(f"test-4nt8-{traj}")
        psel = mk_psel(total_site, rate, rs)
        psel.save(f"{path_src}/traj-{traj}.txt")
    return total_site, path_src, path_dst

@q.timer(is_timer_fork=True)
def run_conversion(total_site, path_dst, path_src):
    geo = q.Geometry(total_site)
    fn_list = glob.glob(f"*.txt", root_dir=path_src)
    fn_list.sort()
    for fn in fn_list:
        name = fn.removesuffix(".txt")
        fn_src = f"{path_src}/{name}.txt"
        fn_dst = f"{path_dst}/{name}.lati"
        psel = q.PointsSelection()
        psel.load(fn_src, geo=geo)
        psel.save(fn_dst)
        if is_test():
            psel_load = q.PointsSelection()
            psel_load.load(fn_dst)
            assert pickle.dumps(psel_load) == pickle.dumps(psel)
            q.json_results_append(f"fn_src={fn_src}")
            q.json_results_append(f"fn_dst={fn_dst}")
            psel_str = f"{psel.total_site} {psel[:].tolist()}"
            q.json_results_append(f"hash(psel)={q.hash_sha256(psel_str)}")
            q.displayln_info(f"psel: {psel_str}")

def is_test():
    return q.get_arg("--src") is None

@q.timer(is_timer_fork=True)
def run():
    total_site_str = q.get_arg("--grid")
    path_src = q.get_arg("--src")
    path_dst = q.get_arg("--dst")
    if is_test():
        assert total_site_str is None
        assert path_src is None
        assert path_dst is None
        q.displayln_info(f"Usage:{usage}")
        q.displayln_info(f"Will now generate test data and run conversion.")
        total_site, path_src, path_dst = gen_test_data()
    else:
        assert isinstance(total_site_str, str)
        assert path_src is not None
        assert path_dst is not None
        total_site = parse_coordinate(total_site_str)
    run_conversion(total_site, path_dst, path_src)

if __name__ == "__main__":
    q.begin_with_mpi()
    run()
    q.check_log_json(__file__, check_eps=1e-10)
    q.end_with_mpi()
    q.displayln_info(f"CHECK: finished successfully.")
