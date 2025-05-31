#!/usr/bin/env python3

import gpt as g
import qlat as q
import qlat_gpt as qg
from qlat_scripts.v1 import (
        get_param,
        set_param,
        mk_inverter,
        )
import subprocess

qg.begin_with_gpt()

job_tag = "test-4nt8"
traj = 1000

q.qremove_all_info("results")
q.qmkdir_info("results")

total_site = q.Coordinate(get_param(job_tag, "total_site"))
geo = q.Geometry(total_site)
q.displayln_info("CHECK: geo.show() =", geo.show())
rs = q.RngState(f"seed-{job_tag}-{traj}")

gf = q.GaugeField(geo)
gf.set_rand(rs.split("gf-init"), 0.05, 2)
gf.unitarize()
gf.show_info()

gt = qg.gauge_fix_coulomb(gf)

n_per_tslice_smear = 8
fsel_smear = q.FieldSelection()
fsel_smear.set_rand(total_site, n_per_tslice_smear, rs)
psel_smear = fsel_smear.to_psel()

n_points = 256
psel = q.PointsSelection()
psel.set_rand(total_site, n_points, rs)

n_per_tslice = 16
fsel = q.FieldSelection()
fsel.set_rand(total_site, n_per_tslice, rs)

fselc = fsel.copy()
fselc.add_psel(psel)

gf_ape = q.gf_spatial_ape_smear(gf, 0.5, 30)
gf_ape = q.mk_left_expanded_gauge_field(gf_ape)

inv_type = 1
inv_acc = 0

inv = mk_inverter(gf, job_tag, inv_type, inv_acc, n_grouped = q.get_num_node())

xg = q.Coordinate(psel_smear[0])

tag = f"smear ; xg={tuple(xg.to_list())} ; type={inv_type} ; accuracy={inv_acc}"

q.displayln_info("CHECK: ", tag)
src = q.mk_point_src(geo, xg)

q.displayln_info(f"CHECK: qnorm(src) = {q.qnorm(src)}")

smear_coef = 0.9375
smear_step = 10
src = q.prop_smear(src, gf_ape, smear_coef, smear_step)

q.displayln_info(f"CHECK: qnorm(src) = {q.qnorm(src):.12E} after smear")

sol = inv * src

q.displayln_info(f"CHECK: qnorm(sol) = {q.qnorm(sol):.6E}")

sol_psel = q.PselProp(psel)
sol_psel @= sol
sol_psel.save(f"results/{tag} ; psnk.lat")

q.displayln_info(f"CHECK: qnorm(sol_psel) = {q.qnorm(sol_psel):.6E}")

sfw = q.open_fields("results/prop-smear", "w", q.Coordinate([ 1, 1, 2, 4, ]))

sol_s = q.SelProp(fselc)
sol_s @= sol

q.displayln_info(f"CHECK: qnorm(sol_s) = {q.qnorm(sol_s):.6E}")

sol_s.save_float_from_double(sfw, f"{tag}")

sfw.close()

sol_gt = gt * sol

sol_ws = sol_gt.glb_sum_tslice()
sol_ws.save(f"results/{tag} ; wsnk.lat")

q.displayln_info(f"CHECK: qnorm(sol_ws) = {q.qnorm(sol_ws):.5E}")

sol_smear_psel = q.PselProp(psel_smear)

sol_smear = q.prop_smear(sol, gf_ape, smear_coef, smear_step)

q.displayln_info(f"CHECK: qnorm(sol_smear) = {q.qnorm(sol_smear):.5E}")

sol_smear_psel = q.PselProp(psel_smear)
sol_smear_psel @= sol_smear
sol_smear_psel.save(f"results/{tag} ; smear-snk.lat")

q.displayln_info(f"CHECK: qnorm(sol_smear_psel) = {q.qnorm(sol_smear_psel):.5E}")

q.check_all_files_crc32_info("results")

q.timer_display()

qg.end_with_gpt()

q.displayln_info(f"CHECK: finished successfully.")
