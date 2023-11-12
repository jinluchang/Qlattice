#!/usr/bin/env python3

import qlat as q
import numpy as np

import sys, os

q.begin_with_mpi()

q.qremove_all_info("results")
q.qmkdir_info("results")
rs = q.RngState("seed")
total_site = q.Coordinate([4, 4, 4, 8])
geo = q.Geometry(total_site, 1)
q.displayln_info("CHECK: geo.show() =", geo.show())

psel = q.PointsSelection([ [0,0,0,0], [0,1,2,0], ])

psel.save("results/psel.txt")

psel = q.PointsSelection()
psel.load("results/psel.txt", geo)

q.displayln_info("CHECK: psel.xg_arr().tolist() =", psel.xg_arr().tolist())

n_per_tslice = 16

fsel = q.FieldSelection()
fsel.set_rand_psel(geo.total_site(), n_per_tslice, rs.split("fsel"))

prob = n_per_tslice * total_site[3] / geo.total_volume()

q.displayln_info("CHECK: fsel info =", fsel.geo().show(), n_per_tslice, f"{prob:.14E}")

total_n_elems = q.glb_sum(fsel.n_elems())
q.displayln_info(f"CHECK: fsel total_n_elems = {total_n_elems}")

fsel.save("results/fsel.field")

fsel = q.FieldSelection()
fsel.load("results/fsel.field")
fsel.save("results/fsel-1.field")

q.displayln_info("CHECK: fsel info =", fsel.geo().show(), n_per_tslice, f"{prob:.14E}")

total_n_elems = q.glb_sum(fsel.n_elems())
q.displayln_info(f"CHECK: fsel total_n_elems = {total_n_elems}")

fsel.add_psel(psel)

q.displayln_info("CHECK: fsel info =", fsel.geo().show(), n_per_tslice, f"{prob:.14E}")

total_n_elems = q.glb_sum(fsel.n_elems())
q.displayln_info(f"CHECK: fsel total_n_elems = {total_n_elems}")

psel = q.PointsSelection()
psel.set_rand(total_site, 32, rs.split("psel"))
q.displayln_info(f"CHECK: psel.n_points() = {psel.n_points()}")

xg_arr_double = np.array(psel.xg_arr(), dtype=float)
sig = q.get_double_sig(xg_arr_double, rs.split("sig"))
q.displayln_info(f"CHECK: psel.xg_arr sig = {sig:.10E}")

psel1 = psel.intersect(fsel)
q.displayln_info(f"CHECK: psel1.n_points() = {psel1.n_points()}")

xg_arr_double = np.array(psel1.xg_arr(), dtype=float)
sig = q.get_double_sig(xg_arr_double, rs.split("sig"))
q.displayln_info(f"CHECK: psel1.xg_arr sig = {sig:.10E}")

fsel_psel = q.FieldSelection()
fsel_psel.set_empty(geo)

if q.get_id_node() == 0:
    q.displayln_info(os.listdir("results"))

q.timer_display()

q.end_with_mpi()

q.displayln_info(f"CHECK: finished successfully.")
