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
geo = q.Geometry(total_site)
q.displayln_info("CHECK: geo.show() =", geo.show())
rs = q.RngState("seed")

gf = q.GaugeField(geo)

gf.set_rand(rs.split("gf-init"), 0.3, 1)

gf.show_info()

q.json_results_append(f"gf", q.get_data_sig_arr(gf, q.RngState(), 3), 1e-12)

mode_fft = 1
fft_f = q.mk_fft(True, is_normalizing=True, is_only_spatial=False, mode_fft=mode_fft)
fft_b = q.mk_fft(False, is_normalizing=True, is_only_spatial=False, mode_fft=mode_fft)


gfm = fft_f * gf

q.json_results_append(f"gfm", q.get_data_sig_arr(gfm, q.RngState(), 3), 1e-12)

gf1 = fft_b * gfm

gf1.show_info()

q.json_results_append(f"gf1", q.get_data_sig_arr(gf1, q.RngState(), 3), 1e-12)

f_factor = q.mk_phase_field(gf.geo, [1, 0, 0, 0,])

gf *= f_factor

gf.show_info()

q.displayln_info(gf.get_elems_xg(q.Coordinate([0, 0, 0, 0,]))[:])
qnorm = q.qnorm(gf.get_elems_xg(q.Coordinate([0, 0, 0, 0,]))[:])
q.json_results_append(f"qnorm", qnorm, 1e-12)

xg1 = rs.split("xg1").c_rand_gen(total_site)
xg2 = rs.split("xg2").c_rand_gen(total_site)

val1 = gf.get_elems_xg(xg1)[:]
q.displayln_info(f"gf xg={xg1} val1={val1}")
sig1 = q.get_data_sig(val1, rs.split("sig1"))
q.displayln_info(f"CHECK: gf xg={xg1}")
q.json_results_append(f"gf sig1", sig1, 1e-12)

val12 = gf.get_elems_xg([ xg1, xg2, ])[:]
q.displayln_info(f"gf xg={xg1, xg2} val12={val12}")
sig12 = q.get_data_sig(val12, rs.split("sig12"))
q.displayln_info(f"CHECK: gf xg={xg1, xg2}")
q.json_results_append(f"gf sig12", sig12, 1e-12)

m = 1
val2 = gf.get_elem_xg(xg1, m)[:]
q.displayln_info(f"gf xg={xg1} val2={val2}")
sig2 = q.get_data_sig(val2, rs.split("sig2"))
q.displayln_info(f"CHECK: gf xg={xg1} m={m}")
q.json_results_append(f"gf sig2", sig2, 1e-12)

val22 = gf.get_elem_xg([ xg1, xg2, ], m)[:]
q.displayln_info(f"gf xg={xg1, xg2} val22={val22}")
sig22 = q.get_data_sig(val22, rs.split("sig22"))
q.displayln_info(f"CHECK: gf xg={xg1, xg2} m={m}")
q.json_results_append(f"gf sig22", sig22, 1e-12)

gf_sum_initial = gf.glb_sum()

gf_sum = gf_sum_initial.copy()[:]

q.displayln_info(gf_sum)
qnorm = q.qnorm(gf_sum)
q.json_results_append(f"qnorm", qnorm, 1e-12)

gf_sum_tslice = gf.glb_sum_tslice()

for t in range(total_site[3]):
    gf_sum -= gf_sum_tslice.get_elems(t)[:]

q.displayln_info(q.qnorm(gf_sum))

assert q.qnorm(gf_sum) <= 1e-16

qnorm = q.qnorm(gf_sum_tslice.to_numpy())
q.json_results_append(f"t_dir=3 (default) qnorm(gf_sum_tslice)", qnorm, 1e-12)

for t_dir in range(4):
    gf_sum_tslice = gf.glb_sum_tslice(t_dir = t_dir)
    n_points = gf_sum_tslice.n_points
    multiplicity = gf_sum_tslice.multiplicity
    psel_list = gf_sum_tslice.psel.xg_arr.tolist()
    q.displayln_info(f"CHECK: t_dir={t_dir} n_points={n_points} multiplicity={multiplicity} psel_list={psel_list}")
    qnorm = q.qnorm(gf_sum_tslice.to_numpy())
    q.json_results_append(f"t_dir={t_dir} qnorm(gf_sum_tslice)", qnorm, 1e-12)
    gf_sum = gf_sum_initial.copy()[:]
    for t in range(total_site[t_dir]):
        gf_sum -= gf_sum_tslice.get_elems(t)[:]
    q.displayln_info(f"t_dir={t_dir} qnorm diff", f"{q.qnorm(gf_sum):.14E}")
    assert q.qnorm(gf_sum) <= 1e-23 * q.qnorm(gf_sum_initial)

f = gf.as_field(q.ElemTypeComplexD)

q.displayln_info(f.get_elems([0, 0, 0, 0,]))
qnorm = q.qnorm(f.get_elems([0, 0, 0, 0,]))
q.json_results_append(f"qnorm", qnorm, 1e-12)

gf1 = q.GaugeField()

gf1.from_field(f)

gf1 -= gf

q.json_results_append(f"diff norm", gf1.qnorm(), 1e-12)

n_per_tslice = 16
fsel = q.FieldSelection()
fsel.set_rand(geo.total_site, n_per_tslice, rs.split("fsel"))

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

q.json_results_append(f"sum_tslice1.qnorm()", sum_tslice1.qnorm(), 1e-12)
q.json_results_append(f"sum_tslice2.qnorm()", sum_tslice2.qnorm(), 1e-12)
q.json_results_append(f"sum_tslice.qnorm()", sum_tslice.qnorm(), 1e-12)

gf_vec = [ q.mk_merged_fields_ms([ (gf, m,), ]) for m in range(4) ]

gf1 = q.GaugeField()

gf1 @= q.mk_merged_fields_ms([ (gf_vec[m], 0) for m in range(4) ])

gf1.show_info()

gf2 = q.GaugeField()

gf2 @= q.mk_merged_fields_ms([ (gf, m) for m in range(4) ])

gf2.show_info()

f = q.FieldRealD(geo, 2)
f.set_rand(q.RngState("seed-f-init"))

sig = q.glb_sum(q.get_data_sig(f[:], q.RngState(f"seed-f-sig-{q.get_id_node()}")))
q.json_results_append(f"f sig", sig, 1e-7)

f1 = f.shift()
sig = q.glb_sum(q.get_data_sig(f1[:], q.RngState(f"seed-f1-sig-{q.get_id_node()}")))
q.json_results_append(f"f1 sig", sig, 1e-7)

f2 = f.shift(q.Coordinate([ 1, 2, 3, 4, ]))
sig = q.glb_sum(q.get_data_sig(f2[:], q.RngState(f"seed-f2-sig-{q.get_id_node()}")))
q.json_results_append(f"f2 sig", sig, 1e-7)

f3 = f.shift(q.Coordinate([ 0, -2, 8, -9, ]), True)
sig = q.glb_sum(q.get_data_sig(f3[:], q.RngState(f"seed-f3-sig-{q.get_id_node()}")))
q.json_results_append(f"f3 sig", sig, 1e-7)

f4 = f.shift(is_reflect=True)
sig = q.glb_sum(q.get_data_sig(f4[:], q.RngState(f"seed-f4-sig-{q.get_id_node()}")))
q.json_results_append(f"f4 sig", sig, 1e-7)

fsel = q.FieldSelection(geo)
fsel[q.RngState("seed-fsel-init").u_rand_arr(fsel[:].shape) < 1/16] = 0
fsel.update()
q.displayln_info(f"CHECK: q.glb_sum(fsel.n_elems) = {q.glb_sum(fsel.n_elems)}")

sf = q.SelectedFieldRealD(fsel, 2)
sf @= f

sig = q.glb_sum(q.get_data_sig(sf[:], q.RngState(f"seed-sf-sig-{q.get_id_node()}")))
q.json_results_append(f"sf sig", sig, 1e-7)

sf1 = sf.shift()
sig = q.glb_sum(q.get_data_sig(sf1[:], q.RngState(f"seed-sf1-sig-{q.get_id_node()}")))
q.json_results_append(f"sf1 sig", sig, 1e-7)

sf1 @= f1
sig = q.glb_sum(q.get_data_sig(sf1[:], q.RngState(f"seed-sf1-sig-{q.get_id_node()}")))
q.json_results_append(f"sf1 sig", sig, 1e-7)

sf2 = sf.shift(q.Coordinate([ 1, 2, 3, 4, ]))
sig = q.glb_sum(q.get_data_sig(sf2[:], q.RngState(f"seed-sf2-sig-{q.get_id_node()}")))
q.json_results_append(f"sf2 sig", sig, 1e-7)

sf2 @= f2
sig = q.glb_sum(q.get_data_sig(sf2[:], q.RngState(f"seed-sf2-sig-{q.get_id_node()}")))
q.json_results_append(f"sf2 sig", sig, 1e-7)

sf3 = sf.shift(q.Coordinate([ 0, -2, 8, -9, ]), True)
sig = q.glb_sum(q.get_data_sig(sf3[:], q.RngState(f"seed-sf3-sig-{q.get_id_node()}")))
q.json_results_append(f"sf3 sig", sig, 1e-7)

sf3 @= f3
sig = q.glb_sum(q.get_data_sig(sf3[:], q.RngState(f"seed-sf3-sig-{q.get_id_node()}")))
q.json_results_append(f"sf3 sig", sig, 1e-7)

sf4 = sf.shift(is_reflect=True)
sig = q.glb_sum(q.get_data_sig(sf4[:], q.RngState(f"seed-sf4-sig-{q.get_id_node()}")))
q.json_results_append(f"sf4 sig", sig, 1e-7)

sf4 @= f4
sig = q.glb_sum(q.get_data_sig(sf4[:], q.RngState(f"seed-sf4-sig-{q.get_id_node()}")))
q.json_results_append(f"sf4 sig", sig, 1e-7)

q.check_log_json(__file__)
q.timer_display()
q.end_with_mpi()
q.displayln_info(f"CHECK: finished successfully.")
