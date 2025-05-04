#!/usr/bin/env python3

import qlat as q
import numpy as np

import sys, os

q.begin_with_mpi()

q.qremove_all_info("results")
q.qmkdir_info("results")
rs = q.RngState("seed")
total_site = q.Coordinate([4, 4, 4, 8])
geo = q.Geometry(total_site)
q.json_results_append(f"geo.show() = {geo.show()}")

psel = q.PointsSelection(total_site, [ [0,0,0,0], [0,1,2,0], ])

psel.save("results/psel.txt")

psel = q.PointsSelection()
psel.load("results/psel.txt", geo)

q.json_results_append(f"psel.xg_arr.tolist() = {psel.xg_arr.tolist()}")

n_per_tslice = 16

fsel = q.FieldSelection()
fsel.set_rand_psel(geo.total_site, n_per_tslice, rs.split("fsel"))

prob = n_per_tslice * total_site[3] / geo.total_volume

total_n_elems = q.glb_sum(fsel.n_elems)
q.json_results_append(f"fsel info = {fsel.geo.show()} total_n_elems = {total_n_elems}")

fsel.save("results/fsel.field")
fsel.load("results/fsel.field")
fsel.save("results/fsel-1.field")

total_n_elems = q.glb_sum(fsel.n_elems)
q.json_results_append(f"fsel info = {fsel.geo.show()} total_n_elems = {total_n_elems}")

fsel.add_psel(psel)

total_n_elems = q.glb_sum(fsel.n_elems)
q.json_results_append(f"fsel info = {fsel.geo.show()} total_n_elems = {total_n_elems}")

psel = q.PointsSelection()
psel.set_rand(total_site, 32, rs.split("psel"))
q.json_results_append(f"psel.n_points = {psel.n_points}")

xg_arr_double = np.array(psel.xg_arr, dtype=float)
sig = q.get_data_sig(xg_arr_double, rs.split("sig"))
q.json_results_append(f"psel.xg_arr sig", sig)

psel1 = psel.intersect(fsel)
q.json_results_append(f"psel1.n_points = {psel1.n_points}")

xg_arr_double = np.array(psel1.xg_arr, dtype=float)
sig = q.get_data_sig(xg_arr_double, rs.split("sig"))
q.json_results_append(f"psel1.xg_arr sig", sig)

fsel2 = q.FieldSelection()
fsel2.set_empty(geo)
fsel2.set_rand(total_site, n_per_tslice, rs.split("fsel2"))

total_n_elems = q.glb_sum(fsel2.n_elems)
q.json_results_append(f"fsel2 info = {fsel2.geo.show()} total_n_elems = {total_n_elems}")

fsel3 = fsel2.copy()

fsel3.add_fsel(fsel)
total_n_elems = q.glb_sum(fsel3.n_elems)
q.json_results_append(f"fsel3 info = {fsel3.geo.show()} total_n_elems = {total_n_elems}")

fsel4 = fsel2.intersect(fsel)
total_n_elems = q.glb_sum(fsel4.n_elems)
q.json_results_append(f"fsel4 info = {fsel4.geo.show()} total_n_elems = {total_n_elems}")

fsel5 = fsel.intersect(fsel2)
total_n_elems = q.glb_sum(fsel5.n_elems)
q.json_results_append(f"fsel5 info = {fsel4.geo.show()} total_n_elems = {total_n_elems}")

if q.get_id_node() == 0:
    q.displayln_info(os.listdir("results"))

q.check_log_json(__file__, check_eps=1e-10)

q.timer_display()

q.end_with_mpi()

q.displayln_info(f"CHECK: finished successfully.")
