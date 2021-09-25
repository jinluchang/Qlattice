#!/usr/bin/env python3

# Need --mpi X.X.X.X runtime option

import qlat_ext as q
import gpt as g
import qlat_gpt as qg
import pprint
import math
import rbc_ukqcd_params as rup
import rbc_ukqcd as ru
import numpy as np

import os

def get_save_path(fn):
    return os.path.join("results", fn)

def get_load_path(fn):
    if fn is None:
        return None
    path_list = [
            "results",
            "/home/frank/application/Public/Muon-GM2-cc/jobs/final-run/all-analysis-data",
            ]
    for path in path_list:
        p = os.path.join(path, fn)
        if q.does_file_exist_sync_node(p):
            return p
    return None

@q.timer
def mk_pion_prop(total_site, pion_mass):
    geo = q.Geometry(total_site, 1)
    grid = qg.mk_grid(geo)
    c = g.complex(grid)
    c[:] = 0
    xg = [ 0, 0, 0, 0, ]
    xl = geo.coordinate_l_from_g(xg)
    c[ [xg,] ] = 1.0
    f = qg.qlat_from_gpt([c])
    f = q.free_scalar_invert_cfield(f, pion_mass)
    for t in range(12):
        for i in range(16):
            q.displayln_info("spatial", t, i, f.get_elems([0, 0, i, t,])[0].real)
    return f

r_scaling_factor = 5.0

@q.timer
def mk_four_point_func_table(total_site, n_dtype):
    info_list = [
            [ "type", n_dtype, ],
            [ "t", math.ceil(total_site[3] / 2), ],
            [ "r", math.ceil(1.0 + r_scaling_factor * math.sqrt(3.0) *
                total_site[0] / 2.0), ],
            [ "em", 4, [ "mm", "tt", "ii", "xx", ], ],
            ]
    ld = q.mk_ld(*info_list)
    return ld

@q.timer
def mk_four_point_func_table_ff(total_site, pion_mass, tag = ""):
    r_pi_list = [ 3.0, ]
    n_dtype = len(r_pi_list)
    ld = mk_four_point_func_table(total_site, n_dtype)
    for dtype, r_pi in enumerate(r_pi_list):
        f = q.mk_pion_four_point_field(total_site, pion_mass, tag, r_pi)
        q.acc_four_point_func_em(ld, f, dtype, r_scaling_factor)
    ld.save(f"results/table-{total_site[0]}nt{total_site[3]}-{pion_mass}-({tag}).lat")
    # q.partial_sum_r_four_point_func_em(ld)

rup.dict_params["48I"]["data_path_em"] = get_load_path("lat-four-point-em/48I-0.00078")
rup.dict_params["48I"]["ainv/gev"] = 1.73

rup.dict_params["64I"]["data_path_em"] = get_load_path("lat-four-point-em/64I-0.000678")
rup.dict_params["64I"]["ainv/gev"] = 2.359

@q.timer
def get_four_point_em(job_tag, traj, dtype):
    path = rup.dict_params[job_tag]["data_path_em"]
    fn = f"{path}/results={traj}/four-point-func-em.lat"
    ld = q.LatData()
    ld.load(fn)
    ld[(0,)] *= 1 / space_volume
    ld[(1,)] *= 1 / space_volume

    if dtype in [ 0, 1, ]:
        total_site = rup.dict_params[job_tag]["total_site"]
        space_volume = total_site[0] * total_site[1] * total_site[2]
    q.partial_sum_r_four_point_func_em(ld)
    r_max = ld.dim_size(2) - 1
    for dtype in [ 0, 1, ]:
        for t in range(ld.dim_size(1)):
            fac = 2 / ld[(1, t, r_max, 1,)]
            if t == 0:
                fac /= 2
            ld[(dtype, t,)] *= fac
    return ld

@q.timer
def run_job(total_site, pion_mass):
    f = mk_pion_prop(total_site, pion_mass)
    fpft = mk_four_point_func_table(total_site, pion_mass)
    fpft = mk_four_point_func_table(total_site, pion_mass, "pole")
    fpft = mk_four_point_func_table(total_site, pion_mass, "0-pole")
    #
    q.clean_cache()
    q.timer_display()

qg.begin_with_gpt()

q.check_time_limit()

total_site = [ 32, 32, 32, 64, ]
pion_mass = 0.135

q.qmkdir_info(f"results")

# run_job(total_site, pion_mass)

q.displayln_info(get_four_point_em("64I", 2810).show())

qg.end_with_gpt()
