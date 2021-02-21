#!/usr/bin/env python3

import qlat as q

import sys, os

q.begin()

rs = q.RngState("seed")

geo = q.Geometry((4, 4, 4, 8), 4)

q.displayln_info("geo.show() =", geo.show())

q.qremove_all_info("results")

q.qmkdir_info("results")
if q.get_id_node() == 0:
    q.displayln_info(os.listdir("results"))

psel = q.PointSelection([(0,0,0,0), (0,1,2,0)])

psel.save("results/psel.txt")

psel = q.PointSelection()
psel.load("results/psel.txt")

q.displayln_info("psel.list() =", psel.list())

fsel = q.FieldSelection(geo.total_site(), 16, rs)

q.displayln_info("fsel info =", fsel.geo().show(), fsel.n_per_tslice(), fsel.prob())

fsel.save("results/fsel.field")

fsel = q.FieldSelection()
fsel.load("results/fsel.field", 16)
fsel.save("results/fsel-1.field")

q.displayln_info("fsel info =", fsel.geo().show(), fsel.n_per_tslice(), fsel.prob())

fsel.add_psel(psel)

q.displayln_info("fsel info =", fsel.geo().show(), fsel.n_per_tslice(), fsel.prob())

if q.get_id_node() == 0:
    q.displayln_info(os.listdir("results"))

q.timer_display()

q.end()
