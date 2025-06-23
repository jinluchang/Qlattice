#!/usr/bin/env python3

import qlat as q

import qlat_gpt as qg
from qlat_scripts.v1 import *

qg.begin_with_gpt()

job_tag = "test-8nt16"
traj = 1000

q.qremove_all_info("results")
q.qmkdir_info("results")

total_site = q.Coordinate(get_param(job_tag, "total_site"))
geo = q.Geometry(total_site)
q.displayln_info("CHECK: geo.show() =", geo.show())
rs = q.RngState(f"seed-{get_job_seed(job_tag)}-{traj}")

gf = q.GaugeField(geo)
gf.set_rand(rs.split("gf-init"), 0.05, 2)
gf.show_info()

inv_type = 1
inv_acc = 0

inv = mk_inverter(gf, job_tag, inv_type, inv_acc, n_grouped = q.get_num_node())

srcs = []

for i in range(2):
    p = q.Prop(geo)
    p.set_rand(rs.split(f"prop-src {i}"))
    q.displayln_info(f"CHECK: prop src {i} qnorm = {p.qnorm():.12E} crc32 = {p.crc32()}")
    srcs.append(p)

sols = inv * srcs

for i in range(2):
    p = sols[i]
    q.displayln_info(f"CHECK: prop sol {i} qnorm = {p.qnorm():.6E}")
    q.displayln_info(f"crc32 = {p.crc32():08X}")

q.timer_display()

qg.end_with_gpt()

q.displayln_info(f"CHECK: finished successfully.")
