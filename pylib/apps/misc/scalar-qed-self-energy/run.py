#!/usr/bin/env python3

# Need --mpi X.X.X.X runtime option

import qlat_ext as q
import gpt as g
import qlat_gpt as qg
import pprint
import math

import os

def get_save_path(fn):
    return os.path.join("results", fn)

def get_load_path(fn):
    if fn is None:
        return None
    path_list = [
            "results",
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
def mk_four_point_func_table(total_site, pion_mass, tag = ""):
    info_list = [
            [ "type", 1, ],
            [ "t", math.ceil(total_site[3] / 2), ],
            [ "r", math.ceil(1.0 + r_scaling_factor * math.sqrt(3.0) *
                total_site[0] / 2.0), ],
            [ "em", 4, [ "mm", "tt", "ii", "xx", ], ],
            ]
    ld = q.mk_ld(*info_list)
    r_pi = 3.0
    f = q.mk_pion_four_point_field(total_site, pion_mass, tag, r_pi)
    q.acc_four_point_func_em(ld, f, 0, r_scaling_factor)
    q.partial_sum_r_four_point_func_em(ld)
    ld.save(f"table({tag}).lat")
    for t in range(12):
        for i in range(16):
            q.displayln_info("spatial", t, i, f.get_elems([0, 0, i, t,])[15].real)
    psf = f.glb_sum_tslice()
    for i in range(32):
        q.displayln_info(i, psf.get_elem(i, 15))

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

run_job(total_site, pion_mass)

qg.end_with_gpt()
