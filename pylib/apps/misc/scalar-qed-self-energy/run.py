#!/usr/bin/env python3

# Need --mpi X.X.X.X runtime option

import qlat_ext as q
import gpt as g
import qlat_gpt as qg
import pprint

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
def set_pion_prop(total_site, pion_mass):
    geo = q.Geometry(total_site, 1)
    # f = q.Field("Complex", geo)
    grid = qg.mk_grid(geo)
    c = g.complex(grid)
    c[:] = 0
    xg = [ 0, 0, 0, 0, ]
    xl = geo.coordinate_l_from_g(xg)
    c[ [xg,] ] = 1.0
    print("before", c[ 0, 0, 0, 0, ], q.get_id_node())
    f = qg.qlat_from_gpt([c])
    f = q.free_scalar_invert_cfield(f, pion_mass)
    c = qg.gpt_from_qlat(f)[0]
    print("after", c[ 0, 0, 0, 0, ], q.get_id_node())
    print("after", c[ 0, 0, 0, 1, ], q.get_id_node())
    return f

@q.timer
def calc_scalar_qed_mom(total_site, pion_mass):
    pass

@q.timer
def run_job(total_site, pion_mass):
    f = set_pion_prop(total_site, pion_mass)
    #
    q.clean_cache()
    q.timer_display()

qg.begin_with_gpt()

q.check_time_limit()

total_site = [ 4, 4, 4, 8, ]
pion_mass = 0.135

run_job(total_site, pion_mass)

qg.end_with_gpt()
