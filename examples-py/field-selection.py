#!/usr/bin/env python3

import qlat as q

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

fsel = q.FieldSelection(geo.total_site(), n_per_tslice, rs)

prob = n_per_tslice * total_site[3] / geo.total_volume()

q.displayln_info("CHECK: fsel info =", fsel.geo().show(), n_per_tslice, f"{prob:.14E}")

fsel.save("results/fsel.field")

fsel = q.FieldSelection()
fsel.load("results/fsel.field")
fsel.save("results/fsel-1.field")

q.displayln_info("CHECK: fsel info =", fsel.geo().show(), n_per_tslice, f"{prob:.14E}")

fsel.add_psel(psel)

q.displayln_info("CHECK: fsel info =", fsel.geo().show(), n_per_tslice, f"{prob:.14E}")

if q.get_id_node() == 0:
    q.displayln_info(os.listdir("results"))

q.timer_display()

q.displayln_info(f"CHECK: finished successfully.")

q.end_with_mpi()
