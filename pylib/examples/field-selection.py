#!/usr/bin/env python3

import qlat as q

import sys, os

q.begin()

q.qremove_all_info("results")
q.qmkdir_info("results")
rs = q.RngState("seed")
total_site = [4, 4, 4, 8]
geo = q.Geometry(total_site, 1)
q.displayln_info("geo.show() =", geo.show())

psel = q.PointSelection([[0,0,0,0], [0,1,2,0]])

psel.save("results/psel.txt")

psel = q.PointSelection()
psel.load("results/psel.txt")

q.displayln_info("psel.to_list() =", psel.to_list())

n_per_tslice = 16

fsel = q.FieldSelection(geo.total_site(), n_per_tslice, rs)

q.displayln_info("fsel info =", fsel.geo().show(), fsel.n_per_tslice(), fsel.prob())

fsel.save("results/fsel.field")

fsel = q.FieldSelection()
fsel.load("results/fsel.field", n_per_tslice)
fsel.save("results/fsel-1.field")

q.displayln_info("fsel info =", fsel.geo().show(), fsel.n_per_tslice(), fsel.prob())

fsel.add_psel(psel)

q.displayln_info("fsel info =", fsel.geo().show(), fsel.n_per_tslice(), fsel.prob())

if q.get_id_node() == 0:
    q.displayln_info(os.listdir("results"))

q.timer_display()

q.end()
