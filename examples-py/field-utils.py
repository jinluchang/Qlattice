#!/usr/bin/env python3

import sys
import qlat as q
import numpy as np

size_node_list = [
    [1, 1, 1, 1],
    [1, 1, 1, 2],
    [1, 1, 1, 4],
    [1, 1, 1, 8],
    [2, 2, 2, 2],
    [2, 2, 2, 4]]

q.begin_with_mpi(size_node_list)

q.qremove_all_info("results")
q.qmkdir_info("results")

total_site = q.Coordinate([ 4, 4, 4, 8, ])
geo = q.Geometry(total_site, 1)
q.displayln_info("CHECK: geo.show() =", geo.show())
rs = q.RngState("seed")

gf = q.GaugeField(geo)

gf.set_rand(rs.split("gf-init"), 0.3, 1)

gf.show_info()

fft_f = q.mk_fft(True, is_normalizing = True, is_only_spatial = False)
fft_b = q.mk_fft(False, is_normalizing = True, is_only_spatial = False)

gfm = fft_f * gf

gf1 = fft_b * gfm

gf1.show_info()

f_factor = q.mk_phase_field(gf.geo(), [1, 0, 0, 0,])

gf *= f_factor

gf.show_info()

q.displayln_info(gf.get_elems_xg(q.Coordinate([0, 0, 0, 0,]))[:])
qnorm = q.qnorm(gf.get_elems_xg(q.Coordinate([0, 0, 0, 0,]))[:])
q.displayln_info(f"CHECK: {qnorm:.14E}")

xg1 = rs.split("xg1").c_rand_gen(total_site)
xg2 = rs.split("xg2").c_rand_gen(total_site)

val1 = gf.get_elems_xg(xg1)[:]
q.displayln_info(f"gf xg={xg1} val1={val1}")
sig1 = q.get_double_sig(val1, rs.split("sig1"))
q.displayln_info(f"CHECK: gf xg={xg1} sig={sig1:.14E}")

val12 = gf.get_elems_xg([ xg1, xg2, ])[:]
q.displayln_info(f"gf xg={xg1, xg2} val12={val12}")
sig12 = q.get_double_sig(val12, rs.split("sig12"))
q.displayln_info(f"CHECK: gf xg={xg1, xg2} sig={sig12:.14E}")

m = 1
val2 = gf.get_elem_xg(xg1, m)[:]
q.displayln_info(f"gf xg={xg1} val2={val2}")
sig2 = q.get_double_sig(val2, rs.split("sig2"))
q.displayln_info(f"CHECK: gf xg={xg1} m={m} sig={sig2:.14E}")

val22 = gf.get_elem_xg([ xg1, xg2, ], m)[:]
q.displayln_info(f"gf xg={xg1, xg2} val22={val22}")
sig22 = q.get_double_sig(val22, rs.split("sig22"))
q.displayln_info(f"CHECK: gf xg={xg1, xg2} m={m} sig={sig22:.14E}")

gf_sum_initial = gf.glb_sum()

gf_sum = gf_sum_initial.copy()[:]

q.displayln_info(gf_sum)
qnorm = q.qnorm(gf_sum)
q.displayln_info(f"CHECK: {qnorm:.14E}")

gf_sum_tslice = gf.glb_sum_tslice()

for t in range(total_site[3]):
    gf_sum -= gf_sum_tslice.get_elems(t)[:]

q.displayln_info(q.qnorm(gf_sum))

assert q.qnorm(gf_sum) <= 1e-16

qnorm = q.qnorm(gf_sum_tslice.to_numpy())
q.displayln_info(f"CHECK: t_dir=3 (default) qnorm(gf_sum_tslice)={qnorm:.14E}")

for t_dir in range(4):
    gf_sum_tslice = gf.glb_sum_tslice(t_dir = t_dir)
    n_points = gf_sum_tslice.n_points()
    multiplicity = gf_sum_tslice.multiplicity()
    psel_list = gf_sum_tslice.psel.xg_arr().tolist()
    q.displayln_info(f"CHECK: t_dir={t_dir} n_points={n_points} multiplicity={multiplicity} psel_list={psel_list}")
    qnorm = q.qnorm(gf_sum_tslice.to_numpy())
    q.displayln_info(f"CHECK: t_dir={t_dir} qnorm(gf_sum_tslice)={qnorm:.14E}")
    gf_sum = gf_sum_initial.copy()[:]
    for t in range(total_site[t_dir]):
        gf_sum -= gf_sum_tslice.get_elems(t)[:]
    q.displayln_info(f"t_dir={t_dir} qnorm diff", f"{q.qnorm(gf_sum):.14E}")
    assert q.qnorm(gf_sum) <= 1e-23 * q.qnorm(gf_sum_initial)

f = gf.as_field(q.ElemTypeComplexD)

q.displayln_info(f.get_elems([0, 0, 0, 0,]))
qnorm = q.qnorm(f.get_elems([0, 0, 0, 0,]))
q.displayln_info(f"CHECK: {qnorm:.14E}")

gf1 = q.GaugeField()

gf1.from_field(f)

gf1 -= gf

q.displayln_info("CHECK: diff norm", gf1.qnorm())

n_per_tslice = 16
fsel = q.FieldSelection()
fsel.set_rand(geo.total_site(), n_per_tslice, rs.split("fsel"))

prop = q.Prop(geo)
prop.set_rand(rs.split("prop-1"))
s_prop = q.SelProp(fsel)
s_prop @= prop
prop.set_zero()
prop @= s_prop

sum_tslice1 = prop.glb_sum_tslice()

sum_tslice2 = s_prop.glb_sum_tslice()

sum_tslice = sum_tslice1.copy()
sum_tslice -= sum_tslice2

q.displayln_info(f"CHECK: {sum_tslice1.qnorm():.14E} {sum_tslice2.qnorm():.14E} {sum_tslice.qnorm():.14E}")

gf_vec = [ q.mk_merged_fields_ms([ (gf, m,), ]) for m in range(4) ]

gf1 = q.GaugeField()

gf1 @= q.mk_merged_fields_ms([ (gf_vec[m], 0) for m in range(4) ])

gf1.show_info()

gf2 = q.GaugeField()

gf2 @= q.mk_merged_fields_ms([ (gf, m) for m in range(4) ])

gf2.show_info()

f = q.FieldRealD(geo, 2)
f.set_rand(q.RngState("seed-f-init"))

sig = q.glb_sum(q.get_double_sig(f[:], q.RngState(f"seed-f-sig-{q.get_id_node()}")))
q.displayln_info(f"CHECK: f sig = {sig:.8F}")

f1 = f.shift()
sig = q.glb_sum(q.get_double_sig(f1[:], q.RngState(f"seed-f1-sig-{q.get_id_node()}")))
q.displayln_info(f"CHECK: f1 sig = {sig:.8F}")

f2 = f.shift(q.Coordinate([ 1, 2, 3, 4, ]))
sig = q.glb_sum(q.get_double_sig(f2[:], q.RngState(f"seed-f2-sig-{q.get_id_node()}")))
q.displayln_info(f"CHECK: f2 sig = {sig:.8F}")

f3 = f.shift(q.Coordinate([ 0, -2, 8, -9, ]), True)
sig = q.glb_sum(q.get_double_sig(f3[:], q.RngState(f"seed-f3-sig-{q.get_id_node()}")))
q.displayln_info(f"CHECK: f3 sig = {sig:.8F}")

f4 = f.shift(is_reflect=True)
sig = q.glb_sum(q.get_double_sig(f4[:], q.RngState(f"seed-f4-sig-{q.get_id_node()}")))
q.displayln_info(f"CHECK: f4 sig = {sig:.8F}")

fsel = q.FieldSelection(geo)
fsel[q.RngState("seed-fsel-init").u_rand_arr(fsel[:].shape) < 1/16] = 0
fsel.update()
q.displayln_info(f"CHECK: q.glb_sum(fsel.n_elems()) = {q.glb_sum(fsel.n_elems())}")

sf = q.SelectedFieldRealD(fsel, 2)
sf @= f

sig = q.glb_sum(q.get_double_sig(sf[:], q.RngState(f"seed-sf-sig-{q.get_id_node()}")))
q.displayln_info(f"CHECK: sf sig = {sig:.8F}")

sf1 = sf.shift()
sig = q.glb_sum(q.get_double_sig(sf1[:], q.RngState(f"seed-sf1-sig-{q.get_id_node()}")))
q.displayln_info(f"CHECK: sf1 sig = {sig:.8F}")

sf1 @= f1
sig = q.glb_sum(q.get_double_sig(sf1[:], q.RngState(f"seed-sf1-sig-{q.get_id_node()}")))
q.displayln_info(f"CHECK: sf1 sig = {sig:.8F}")

sf2 = sf.shift(q.Coordinate([ 1, 2, 3, 4, ]))
sig = q.glb_sum(q.get_double_sig(sf2[:], q.RngState(f"seed-sf2-sig-{q.get_id_node()}")))
q.displayln_info(f"CHECK: sf2 sig = {sig:.8F}")

sf2 @= f2
sig = q.glb_sum(q.get_double_sig(sf2[:], q.RngState(f"seed-sf2-sig-{q.get_id_node()}")))
q.displayln_info(f"CHECK: sf2 sig = {sig:.8F}")

sf3 = sf.shift(q.Coordinate([ 0, -2, 8, -9, ]), True)
sig = q.glb_sum(q.get_double_sig(sf3[:], q.RngState(f"seed-sf3-sig-{q.get_id_node()}")))
q.displayln_info(f"CHECK: sf3 sig = {sig:.8F}")

sf3 @= f3
sig = q.glb_sum(q.get_double_sig(sf3[:], q.RngState(f"seed-sf3-sig-{q.get_id_node()}")))
q.displayln_info(f"CHECK: sf3 sig = {sig:.8F}")

sf4 = sf.shift(is_reflect=True)
sig = q.glb_sum(q.get_double_sig(sf4[:], q.RngState(f"seed-sf4-sig-{q.get_id_node()}")))
q.displayln_info(f"CHECK: sf4 sig = {sig:.8F}")

sf4 @= f4
sig = q.glb_sum(q.get_double_sig(sf4[:], q.RngState(f"seed-sf4-sig-{q.get_id_node()}")))
q.displayln_info(f"CHECK: sf4 sig = {sig:.8F}")

q.timer_display()

q.end_with_mpi()

q.displayln_info(f"CHECK: finished successfully.")
